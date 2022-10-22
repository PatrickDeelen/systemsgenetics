package org.molgenis.genotype.vcf;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.*;

import net.sf.samtools.util.BlockCompressedInputStream;

import org.apache.commons.lang3.StringUtils;
import org.molgenis.genotype.AbstractRandomAccessGenotypeData;
import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.GenotypeDataException;
import org.molgenis.genotype.Sample;
import org.molgenis.genotype.Sequence;
import org.molgenis.genotype.SimpleSequence;
import org.molgenis.genotype.annotation.Annotation;
import org.molgenis.genotype.annotation.SampleAnnotation;
import org.molgenis.genotype.annotation.VcfAnnotation;
import org.molgenis.genotype.bgen.BgenGenotypeData;
import org.molgenis.genotype.tabix.TabixFileNotFoundException;
import org.molgenis.genotype.tabix.TabixIndex;
import org.molgenis.genotype.tabix.TabixIndex.TabixIterator;
import org.molgenis.genotype.util.CalledDosageConvertor;
import org.molgenis.genotype.util.FixedSizeIterable;
import org.molgenis.genotype.util.ProbabilitiesConvertor;
import org.molgenis.genotype.variant.GeneticVariant;
import org.molgenis.genotype.variant.GeneticVariantMeta;
import org.molgenis.genotype.variant.GenotypeRecord;
import org.molgenis.genotype.variant.ReadOnlyGeneticVariant;
import org.molgenis.genotype.variant.sampleProvider.SampleVariantUniqueIdProvider;
import org.molgenis.genotype.variant.sampleProvider.SampleVariantsProvider;
import org.molgenis.genotype.vcf.VcfGenotypeField.VcfGenotypeFormat;
import org.molgenis.genotype.vcf.VcfGenotypeField.VcfGenotypeFormatSupplier;
import org.molgenis.vcf.VcfInfo;
import org.molgenis.vcf.VcfReader;
import org.molgenis.vcf.VcfRecord;
import org.molgenis.vcf.VcfSample;
import org.molgenis.vcf.meta.VcfMeta;
import org.molgenis.vcf.meta.VcfMetaContig;
import org.molgenis.vcf.meta.VcfMetaInfo;

import com.google.common.base.Function;
import com.google.common.collect.Iterators;
import org.apache.commons.io.IOUtils;
import org.molgenis.genotype.Allele;
import org.molgenis.genotype.variant.sampleProvider.CachedSampleVariantProvider;

public class VcfGenotypeData extends AbstractRandomAccessGenotypeData implements SampleVariantsProvider {

    private static final org.apache.log4j.Logger LOG = org.apache.log4j.Logger.getLogger(VcfGenotypeData.class);
    private final File bzipVcfFile;
    private final TabixIndex tabixIndex;
    private final int sampleVariantProviderUniqueId;
    private final SampleVariantsProvider variantProvider;
    private final VcfMeta vcfMeta;
    private transient Map<String, Annotation> cachedSampleAnnotationsMap;
    private transient GeneticVariant cachedGeneticVariant;
    private transient VcfRecord cachedVcfRecord;
    private static int totalRandomAccessRequest = 0;
    private static int currentlyOpenFileHandlers = 0;
    private static int closedFileHandlers = 0;
    private final double minimumPosteriorProbabilityToCall;
    private VcfGenotypeFormatSupplier genotypeFormatSupplier;
    private final LinkedHashSet<VcfGenotypeFormat> genotypeProbabilitiesFieldPrecedence;
    private final LinkedHashSet<VcfGenotypeFormat> genotypeCallFieldPrecedence;
    private final LinkedHashSet<VcfGenotypeFormat> genotypeDosageFieldPrecedence;


    /**
     * VCF genotype reader
     *
     * @param bzipVcfFile
     * @throws IOException
     * @throws FileNotFoundException
     */
    //public VcfGenotypeData(File bzipVcfFile, double minimumPosteriorProbabilityToCall) throws FileNotFoundException, IOException {
    //	this(bzipVcfFile, 100, minimumPosteriorProbabilityToCall);
    //}
    public VcfGenotypeData(File bzipVcfFile, int cacheSize, double minimumPosteriorProbabilityToCall) throws FileNotFoundException, IOException {
        this(bzipVcfFile, new File(bzipVcfFile.getAbsolutePath() + ".tbi"), cacheSize, minimumPosteriorProbabilityToCall);
    }

    /**
     * VCF genotype reader with default cache of 100
     *
     * @param bzipVcfFile
     * @param tabixIndexFile
     * @throws IOException
     * @throws FileNotFoundException
     */
    public VcfGenotypeData(File bzipVcfFile, File tabixIndexFile, double minimumPosteriorProbabilityToCall) throws FileNotFoundException, IOException {
        this(bzipVcfFile, tabixIndexFile, 100, minimumPosteriorProbabilityToCall);
    }

    public VcfGenotypeData(File bzipVcfFile, File tabixIndexFile, int cacheSize, double minimumPosteriorProbabilityToCall) throws FileNotFoundException,
            IOException {

        if (!bzipVcfFile.isFile()) {
            throw new FileNotFoundException("VCF file not found at " + bzipVcfFile.getAbsolutePath());
        }

        if (!bzipVcfFile.canRead()) {
            throw new IOException("Cannot access VCF file at: " + bzipVcfFile.getAbsolutePath());
        }

        if (!tabixIndexFile.isFile()) {
            throw new TabixFileNotFoundException(tabixIndexFile.getAbsolutePath(), "VCF tabix file not found at " + tabixIndexFile.getAbsolutePath());
        }

        if (!tabixIndexFile.canRead()) {
            throw new IOException("Cannot read tabix file for VCF at: " + tabixIndexFile.getAbsolutePath());
        }

        if (minimumPosteriorProbabilityToCall < 0 || minimumPosteriorProbabilityToCall > 1) {
            throw new GenotypeDataException("Min posterior probability to call must be >0 and <=1 not:" + minimumPosteriorProbabilityToCall);
        }

        this.bzipVcfFile = bzipVcfFile;
        this.tabixIndex = new TabixIndex(tabixIndexFile, bzipVcfFile, null);
        this.minimumPosteriorProbabilityToCall = minimumPosteriorProbabilityToCall;

        try (VcfReader vcfReader = new VcfReader(new BlockCompressedInputStream(bzipVcfFile))) {
            this.vcfMeta = vcfReader.getVcfMeta();
        }

        if (cacheSize <= 0) {
            variantProvider = this;
        } else {
            variantProvider = new CachedSampleVariantProvider(this, cacheSize);
        }

        sampleVariantProviderUniqueId = SampleVariantUniqueIdProvider.getNextUniqueId();

        genotypeProbabilitiesFieldPrecedence =
                new LinkedHashSet<>(Arrays.asList(VcfGenotypeFormat.GP, VcfGenotypeFormat.GT, VcfGenotypeFormat.DS));
        genotypeCallFieldPrecedence =
                new LinkedHashSet<>(Arrays.asList(VcfGenotypeFormat.GT, VcfGenotypeFormat.GP, VcfGenotypeFormat.DS));
        genotypeDosageFieldPrecedence =
                new LinkedHashSet<>(Arrays.asList(VcfGenotypeFormat.DS, VcfGenotypeFormat.GP, VcfGenotypeFormat.GT));

        genotypeFormatSupplier = new VcfGenotypeFormatSupplier();
    }

    @Override
    public Iterator<GeneticVariant> iterator() {
        final BlockCompressedInputStream inputStream;
        try {
            inputStream = new BlockCompressedInputStream(bzipVcfFile);
        } catch (IOException e) {
            throw new GenotypeDataException(e);
        }

        final VcfReader vcfReader = new VcfReader(inputStream);
        Iterator<GeneticVariant> iterator = new Iterator<GeneticVariant>() {
            private final Iterator<VcfRecord> it = vcfReader.iterator();

            @Override
            public boolean hasNext() {
                boolean hasNext = it.hasNext();
                if (!hasNext) {
                    // close vcf reader
                    try {
                        vcfReader.close();
                    } catch (IOException e) {
                        throw new RuntimeException(e);
                    }
                }
                return hasNext;
            }

            @Override
            public GeneticVariant next() {
                VcfRecord vcfRecord = it.next();
                return toGeneticVariant(vcfRecord);
            }

            @Override
            public void remove() {
                it.remove();
            }
        };
        iterator.hasNext();//needed to properly initiate
        return iterator;
    }

    @Override
    public List<Alleles> getSampleVariants(final GeneticVariant variant) {
        // get vcf record for variant
        VcfRecord vcfRecord = getVcfRecord(variant);
        int nrSamples = vcfRecord.getNrSamples();
        if (nrSamples == 0) {
            return Collections.emptyList();
        }

        VcfGenotypeFormat genotypeFormat = genotypeFormatSupplier.getVcfGenotypeFormat(vcfRecord,
                genotypeCallFieldPrecedence);

        if (VcfGenotypeFormat.GT.equals(genotypeFormat)) {
            return getCalledAlleles(variant, vcfRecord);

        } else if (VcfGenotypeFormat.GP.equals(genotypeFormat)) {

            return ProbabilitiesConvertor.convertProbabilitiesToAlleles(
                    getSampleProbilities(variant),
                    variant.getVariantAlleles(),
                    minimumPosteriorProbabilityToCall);

        } else if (VcfGenotypeFormat.DS.equals(genotypeFormat)) {

            return CalledDosageConvertor.convertDosageToAlleles(getSampleDosage(variant), variant.getVariantAlleles());

        } else {

            ArrayList<Alleles> sampleAlleles = new ArrayList<Alleles>(nrSamples);
            for (int i = 0; i < nrSamples; ++i) {
                sampleAlleles.add(Alleles.BI_ALLELIC_MISSING);
            }
            return sampleAlleles;
        }

    }

    private List<Alleles> getCalledAlleles(GeneticVariant variant, VcfRecord vcfRecord) {
        // convert vcf sample alleles to Alleles§
        List<Alleles> alleles = new ArrayList<>(vcfRecord.getNrSamples());

        try {
            for (VcfSample vcfSample : vcfRecord.getSamples()) {
                List<Allele> vcfAlleles = vcfSample.getAlleles();
                alleles.add(Alleles.createAlleles(vcfAlleles));
            }
        } catch (NumberFormatException ex) {
            throw new GenotypeDataException("Error parsing variant: " + variant.getPrimaryVariantId() + " at " + variant.getSequenceName() + ":" + variant.getStartPos(), ex);
        }
        return alleles;
    }

    @Override
    public Map<String, Annotation> getVariantAnnotationsMap() {
        if (cachedSampleAnnotationsMap == null) {
            cachedSampleAnnotationsMap = new LinkedHashMap<String, Annotation>();

            for (VcfMetaInfo info : vcfMeta.getInfoMeta()) {
                cachedSampleAnnotationsMap.put(info.getId(), VcfAnnotation.fromVcfInfo(info));
            }
        }

        return cachedSampleAnnotationsMap;
    }

    @Override
    public List<Sequence> getSequences() {
        List<String> seqNames = getSeqNames();

        // get sequence length by sequence name
        Map<String, Integer> sequenceLengthById = new HashMap<String, Integer>();
        for (VcfMetaContig contig : vcfMeta.getContigMeta()) {
            sequenceLengthById.put(contig.getId(), contig.getLength());
        }

        List<Sequence> sequences = new ArrayList<Sequence>(seqNames.size());
        for (String seqName : seqNames) {
            sequences.add(new SimpleSequence(seqName, sequenceLengthById.get(seqName), this));
        }

        return sequences;
    }

    @Override
    public List<Sample> getSamples() {
        List<Sample> samples = new ArrayList<Sample>();
        for (String sampleName : vcfMeta.getSampleNames()) {
            samples.add(new Sample(sampleName, null, null));
        }
        return samples;
    }

    @Override
    public int cacheSize() {
        return 0;
    }

    @Override
    public List<Boolean> getSamplePhasing(GeneticVariant variant) {
        VcfRecord vcfRecord = getVcfRecord(variant);

        final int nrSamples = vcfRecord.getNrSamples();
        if (nrSamples == 0) {
            return Collections.emptyList();
        }

        List<Boolean> phasing = new ArrayList<>(nrSamples);
        for (VcfSample vcfSample : vcfRecord.getSamples()) {

            List<Boolean> genotypePhasings = vcfSample.getPhasings();

            if (genotypePhasings == null || genotypePhasings.isEmpty()) {
                phasing.add(Boolean.FALSE);
            } else if (genotypePhasings.size() == 1) {
                phasing.add(genotypePhasings.get(0));
            } else if (genotypePhasings.contains(Boolean.FALSE)) {
                phasing.add(Boolean.FALSE);
            } else {
                phasing.add(Boolean.TRUE);
            }

        }
        return phasing;
    }

    @Override
    public boolean arePhasedProbabilitiesPresent(GeneticVariant variant) {
        VcfRecord vcfRecord = getVcfRecord(variant);

        final int nrSamples = vcfRecord.getNrSamples();
        if (nrSamples == 0) {
            return false;
        }

        LinkedHashSet<VcfGenotypeFormat> haplotypeProbabilitiesFields = getVcfHaplotypeFormats(variant);

        // If the requested format is set and present for this variant base decision on this format
        VcfGenotypeFormat genotypeFormat = genotypeFormatSupplier.getVcfGenotypeFormat(
                vcfRecord, haplotypeProbabilitiesFields);

        return (genotypeFormat != null);
    }

    private LinkedHashSet<VcfGenotypeFormat> getVcfHaplotypeFormats(GeneticVariant variant) {
        LinkedHashSet<VcfGenotypeFormat> haplotypeProbabilitiesFields =
                new LinkedHashSet<>(Arrays.asList(VcfGenotypeFormat.HP, VcfGenotypeFormat.ADS));

        if (variant.hasPhasedGenotypes()) {
            haplotypeProbabilitiesFields.add(VcfGenotypeFormat.GT);
        }
        return haplotypeProbabilitiesFields;
    }

    @Override
    public int getSampleVariantProviderUniqueId() {
        return sampleVariantProviderUniqueId;
    }

    @Override
    public Map<String, SampleAnnotation> getSampleAnnotationsMap() {
        return Collections.emptyMap();
    }

    @Override
    public byte[] getSampleCalledDosage(GeneticVariant variant) {
        return CalledDosageConvertor.convertCalledAllelesToCalledDosage(getSampleVariants(variant),
                variant.getVariantAlleles(), variant.getRefAllele());
    }

    @Override
    public float[] getSampleDosage(GeneticVariant variant) {
        VcfRecord vcfRecord = getVcfRecord(variant);

        final int nrSamples = vcfRecord.getNrSamples();
        if (nrSamples == 0) {
            return new float[0];
        }

        float[] dosages;

        VcfGenotypeFormat genotypeFormat = genotypeFormatSupplier.getVcfGenotypeFormat(vcfRecord,
                genotypeDosageFieldPrecedence);

        int idx = getIndexOfGenotypeFormat(vcfRecord, genotypeFormat);

        if (VcfGenotypeFormat.DS.equals(genotypeFormat)) {
            // retrieve sample dosage from sample info
            dosages = new float[nrSamples];
            int i = 0;
            for (VcfSample vcfSample : vcfRecord.getSamples()) {
                String dosage = vcfSample.getData(idx);
                if (dosage == null) {
                    //throw new GenotypeDataException("Missing DS format value for sample [" + vcfMeta.getSampleName(i) + "] at variant [" + variant.getPrimaryVariantId() + "]");
                    dosages[i++] = -1;
                } else {
                    try {
                        //Math abs to prevent -0 due to rounding
                        dosages[i++] = Math.abs((Float.parseFloat(dosage) - 2) * -1);
                    } catch (NumberFormatException e) {
                        throw new GenotypeDataException("Error in sample dosage (DS) value for sample [" + vcfMeta.getSampleName(i) + "], found value: " + dosage);
                    }
                }

            }
        } else if (VcfGenotypeFormat.GP.equals(genotypeFormat)) {
            dosages = ProbabilitiesConvertor.convertProbabilitiesToDosage(getSampleProbilities(variant), minimumPosteriorProbabilityToCall);

        } else if (VcfGenotypeFormat.GT.equals(genotypeFormat)) {

            // calculate sample dosage from called alleles
            dosages = CalledDosageConvertor.convertCalledAllelesToDosage(getSampleVariants(variant),
                    variant.getVariantAlleles(), variant.getRefAllele());

        } else {
            dosages = new float[nrSamples];
            for (int i = 0; i < nrSamples; ++i) {
                dosages[i] = -1;
            }
        }
        return dosages;
    }

    private int getIndexOfGenotypeFormat(VcfRecord vcfRecord, VcfGenotypeFormat genotypeFormat) {
        int formatIndex = -1;
        if (genotypeFormat != null) {
            String genotypeFormatIdentifier = genotypeFormatSupplier.getGenotypeFormatIdentifier(genotypeFormat);
            formatIndex = vcfRecord.getFormatIndex(genotypeFormatIdentifier);
        }
        return formatIndex;
    }

    @Override
    public void close() throws IOException {
        // noop
    }

    @Override
    public boolean isOnlyContaingSaveProbabilityGenotypes() {
        return false;
    }

    @Override
    public float[][] getSampleProbilities(GeneticVariant variant) {
        VcfRecord vcfRecord = getVcfRecord(variant);

        final int nrSamples = vcfRecord.getNrSamples();
        if (nrSamples == 0) {
            return new float[0][0];
        }

        int numberOfAlleles = variant.getAlleleCount();

        float[][] probs;

        // Select the preferred genotype field
        VcfGenotypeFormat genotypeFormat = genotypeFormatSupplier.getVcfGenotypeFormat(
                vcfRecord, genotypeProbabilitiesFieldPrecedence);

        int idx = getIndexOfGenotypeFormat(vcfRecord, genotypeFormat);

        if (VcfGenotypeFormat.GP.equals(genotypeFormat)) {
            // retrieve sample probabilities from sample info
            probs = new float[nrSamples][3];
            int i = 0;
            for (VcfSample vcfSample : vcfRecord.getSamples()) {
                String probabilitiesStr = vcfSample.getData(idx);
                if (probabilitiesStr == null) {
                    //throw new GenotypeDataException("Missing GP format value for sample [" + vcfMeta.getSampleName(i) + "]");
                    probs[i] = new float[]{0, 0, 0};
                } else {
                    if (probabilitiesStr.matches(".*,+\\.,+.*")) {
//                        System.out.println(probabilitiesStr);
                        probabilitiesStr = probabilitiesStr.replaceAll("\\.", "0");
//                        System.out.println(probabilitiesStr);
                    }
                    String[] probabilities = StringUtils.split(probabilitiesStr, ',');
                    if (probabilities.length != 3 || numberOfAlleles != 2) {
//                        throw new GenotypeDataException("Error in sample prob (GP) value for sample [" + vcfMeta.getSampleName(i) + "], found value: " + probabilitiesStr);
                        probs[i] = new float[3];
                    } else {
                        for (int j = 0; j < 3; ++j) {
                            try {
                                probs[i][j] = Float.parseFloat(probabilities[j]);
                            } catch (NumberFormatException e) {
                                throw new GenotypeDataException("Error in sample prob (GP) value for sample [" + vcfMeta.getSampleName(i) + "], found value: " + probabilitiesStr);
                            }
                        }
                    }
                }
                ++i;
            }

        } else if (VcfGenotypeFormat.GT.equals(genotypeFormat)) {

            probs = ProbabilitiesConvertor.convertCalledAllelesToProbability(
                    getSampleVariants(variant), variant.getVariantAlleles());

        } else if (VcfGenotypeFormat.DS.equals(genotypeFormat)) {

            // calculate sample probabilities from sample dosage
            probs = ProbabilitiesConvertor.convertDosageToProbabilityHeuristic(getSampleDosage(variant));

        } else {
            probs = new float[nrSamples][3];
        }
        return probs;
    }

    @Override
    public double[][] getSampleProbabilitiesComplex(GeneticVariant variant) {
        VcfRecord vcfRecord = getVcfRecord(variant);

        final int nrSamples = vcfRecord.getNrSamples();
        if (nrSamples == 0) {
            return new double[0][0];
        }

        int numberOfAlleles = variant.getAlleleCount();

        double[][] probs;

        // Select the preferred genotype field
        VcfGenotypeFormat genotypeFormat = genotypeFormatSupplier.getVcfGenotypeFormat(
                vcfRecord, genotypeProbabilitiesFieldPrecedence);

        int idx = getIndexOfGenotypeFormat(vcfRecord, genotypeFormat);

        if (VcfGenotypeFormat.GP.equals(genotypeFormat)) {

            probs = readGenotypeProbabilitiesComplex(idx, vcfRecord, numberOfAlleles,
                    getCalledAlleles(variant, vcfRecord));

        } else if (VcfGenotypeFormat.GT.equals(genotypeFormat)) {

            probs = ProbabilitiesConvertor.convertCalledAllelesToComplexProbabilities(
                    getCalledAlleles(variant, vcfRecord),
                    variant.getVariantAlleles());

        } else if (VcfGenotypeFormat.DS.equals(genotypeFormat)) {

            // calculate sample probabilities from sample dosage
            probs = ProbabilitiesConvertor.convertProbabilitiesToComplexProbabilities(
                    ProbabilitiesConvertor.convertDosageToProbabilityHeuristic(
                            getSampleDosage(variant)));

        } else {

            probs = new double[nrSamples][3];
        }
        return probs;
    }

    private double[][] readGenotypeProbabilitiesComplex(int formatIndex, VcfRecord vcfRecord,
                                                        int numberOfAlleles, List<Alleles> sampleVariants) {
        int nrSamples = vcfRecord.getNrSamples();
        double[][] probs;// retrieve sample probabilities from sample info
        probs = new double[nrSamples][];
        int sampleIndex = 0;
        for (VcfSample vcfSample : vcfRecord.getSamples()) {
            String probabilitiesStr = vcfSample.getData(formatIndex);
            if (probabilitiesStr == null) {
                //throw new GenotypeDataException("Missing GP format value for sample [" + vcfMeta.getSampleName(sampleIndex) + "]");
                probs[sampleIndex] = new double[0];
            } else {

                String[] probabilities = StringUtils.split(probabilitiesStr, ',');
                // The number of probabilities should match with the number of probabilities expected with
                // the given number of alleles for this variant and the number of called alleles.
                // Check this if there are called sample variants available.
                if (sampleVariants != null) {
                    // Calculate the number of expected probabilities given the
                    // ploidy and the number of alleles for this variant.
                    int calledAlleleCount = sampleVariants.get(sampleIndex).getAlleleCount();
                    int numberOfExpectedProbabilities = BgenGenotypeData
                            .numberOfProbabilitiesForPloidyAlleleCountCombination(calledAlleleCount,
                                    numberOfAlleles - 1);

                    // Test if this corresponds to the actual number of probabilities.
                    if (probabilities.length != numberOfExpectedProbabilities) {
                        throw new GenotypeDataException(String.format(
                                "Error in sample prob (GP) value for sample [%s], found %d values, " +
                                        "while %d were expected based on %d called alleles and %d variant alleles",
                                vcfMeta.getSampleName(sampleIndex), probabilities.length,
                                numberOfExpectedProbabilities, calledAlleleCount, numberOfAlleles));
                    }
                }

                // Initialize a double array for the probabilities for this sample.
                probs[sampleIndex] = new double[probabilities.length];

                // Loop through the probabilities, assigning the present probabilities.
                for (int i = 0; i < probabilities.length; ++i) {

                    String probabilityAsString = probabilities[i];
                    // Assign the present probability, or 0 if this probability is missing.
                    if (probabilityAsString.equals(".")) {
                        probs[sampleIndex][i] = 0;
                    } else {
                        try {
                            // Parse this probability to a double if it is not missing.
                            probs[sampleIndex][i] = Double.parseDouble(probabilityAsString);

                        } catch (NumberFormatException e) {

                            // Throw an exception if this probability is neither missing (denoted by a dot,
                            // nor parsable to a double.
                            throw new GenotypeDataException(String.format(
                                    "Error in genotype probabilities (GP) value for sample [%s], " +
                                            "found value '%s', while a double was expected",
                                    vcfMeta.getSampleName(sampleIndex),
                                    probabilityAsString));
                        }
                    }
                }
            }
            ++sampleIndex;
        }
        return probs;
    }

    @Override
    public double[][][] getSampleProbabilitiesPhased(GeneticVariant variant) {
        VcfRecord vcfRecord = getVcfRecord(variant);

        final int nrSamples = vcfRecord.getNrSamples();
        if (nrSamples == 0) {
            return new double[0][0][0];
        }

        int numberOfAlleles = variant.getAlleleCount();

        double[][][] probabilities;

        VcfGenotypeFormat genotypeFormat = genotypeFormatSupplier.getVcfGenotypeFormat(vcfRecord,
                getVcfHaplotypeFormats(variant));

        int idx = getIndexOfGenotypeFormat(vcfRecord, genotypeFormat);

        if (VcfGenotypeFormat.HP.equals(genotypeFormat)) {
            // retrieve sample probabilities from sample info
            probabilities = new double[nrSamples][][];
            int sampleIndex = 0;
            for (VcfSample vcfSample : vcfRecord.getSamples()) {
                String probabilitiesStr = vcfSample.getData(idx);
                if (probabilitiesStr == null) {
                    //throw new GenotypeDataException("Missing GP format value for sample [" + vcfMeta.getSampleName(sampleIndex) + "]");
                    probabilities[sampleIndex] = new double[2][2];
                } else {
                    // There should be the same number of values as alleles. These should be split by a comma (",")
                    String[] probabilitiesAsStrings = StringUtils.split(probabilitiesStr, ',');
                    // Check if the expected number of values corresponds to the actual number of values,
                    // and throw an exception if this is not the case.
                    if ((probabilitiesAsStrings.length % numberOfAlleles) != 0) {
                        throw new GenotypeDataException(String.format(
                                "Error in haplotype probabilities (HP) value for sample [%s], found %d value(s) (%s), " +
                                        "while a multiple of %d were expected based on total allele count",
                                vcfMeta.getSampleName(sampleIndex), probabilitiesAsStrings.length,
                                probabilitiesStr, numberOfAlleles));
                    }
                    // Should have a multiple of the number of alleles
                    // For every haplotype, allele, get the probabilities and assign these to the probabilities
                    int numberOfHaplotypes = probabilitiesAsStrings.length / numberOfAlleles;
                    probabilities[sampleIndex] = new double[numberOfHaplotypes][numberOfAlleles];
                    // Loop through haplotypes
                    for (int haplotypeIndex = 0; haplotypeIndex < numberOfHaplotypes; haplotypeIndex++) {
                        // Loop through alleles
                        for (int alleleIndex = 0; alleleIndex < numberOfAlleles; alleleIndex++) {
                            // Get the probability as a string corresponding to this haplotype and allele.
                            String probabilityAsString = probabilitiesAsStrings
                                    [haplotypeIndex * numberOfAlleles + alleleIndex];
                            if (probabilityAsString.equals(".")) {
                                probabilities[sampleIndex][haplotypeIndex][alleleIndex] = 0;
                            } else {
                                try {
                                    probabilities[sampleIndex][haplotypeIndex][alleleIndex] =
                                            Double.parseDouble(probabilityAsString);
                                } catch (NumberFormatException e) {
                                    throw new GenotypeDataException(String.format(
                                            "Error in haplotype probabilities (HP) value for sample [%s], " +
                                                    "found value '%s', while a double was expected",
                                            vcfMeta.getSampleName(sampleIndex),
                                            probabilityAsString));
                                }
                            }
                        }
                    }
                }
                ++sampleIndex;
            }
        } else if (VcfGenotypeFormat.ADS.equals(genotypeFormat)) {

            // check if the number of alleles is 2
            if (variant.getAlleleCount() != 2) {
                throw new GenotypeDataException(String.format(
                        "Error in per-haplotype allele dosage (ADS) for variant %s" +
                                "found %d alleles while 2 were expected",
                        variant.getPrimaryVariantId(), numberOfAlleles));
            }

            // retrieve sample probabilities from sample info
            double[][] haplotypeDosages = new double[nrSamples][];
            int sampleIndex = 0;
            for (VcfSample vcfSample : vcfRecord.getSamples()) {
                String probabilitiesStr = vcfSample.getData(idx);
                if (probabilitiesStr == null) {
                    //throw new GenotypeDataException("Missing GP format value for sample [" + vcfMeta.getSampleName(sampleIndex) + "]");
                    haplotypeDosages[sampleIndex] = new double[2];
                } else {
                    // There should be the same number of values as haplotypes. These should be split by a comma (",")
                    String[] haplotypesAsStrings = StringUtils.split(probabilitiesStr, ',');
                    // Check if the expected number of values corresponds to the actual number of values,
                    // and throw an exception if this is not the case.

                    int numberOfHaplotypes = 2;
                    if (haplotypesAsStrings.length != numberOfHaplotypes) {
                        throw new GenotypeDataException(String.format(
                                "Error in per-haplotype allele dosage (ADS) value for sample [%s], " +
                                        "found %d value(s) (%s), while 2 were expected",
                                vcfMeta.getSampleName(sampleIndex), haplotypesAsStrings.length,
                                probabilitiesStr));
                    }
                    // Should have a multiple of the number of alleles
                    // For every haplotype, allele, get the probabilities and assign these to the probabilities
                    haplotypeDosages[sampleIndex] = new double[numberOfHaplotypes];
                    // Loop through haplotypes
                    for (int haplotypeIndex = 0; haplotypeIndex < numberOfHaplotypes; haplotypeIndex++) {
                        String haplotypeDosageAsString = haplotypesAsStrings[haplotypeIndex];
                        if (haplotypeDosageAsString.equals(".")) {
                            haplotypeDosages[sampleIndex][haplotypeIndex] = 0;
                        } else {
                            try {
                                haplotypeDosages[sampleIndex][haplotypeIndex] = Double.parseDouble(haplotypeDosageAsString);
                            } catch (NumberFormatException e) {
                                throw new GenotypeDataException(String.format(
                                        "Error in per-haplotype allele dosage (ADS) value for sample [%s], " +
                                                "found value '%s', while a double was expected",
                                        vcfMeta.getSampleName(sampleIndex), haplotypeDosageAsString));
                            }
                        }
                    }
                }
                ++sampleIndex;
            }
            probabilities = ProbabilitiesConvertor.haplotypeDosagesToHaplotypeProbabilities(haplotypeDosages);
        } else if (VcfGenotypeFormat.GT.equals(genotypeFormat)) {
            // calculate sample probabilities from sample dosage
            probabilities = ProbabilitiesConvertor.convertCalledAllelesToPhasedProbabilities(getSampleVariants(variant), variant.getVariantAlleles());
        } else {
            throw new GenotypeDataException("Phased data not available");
        }
        return probabilities;
    }

    @Override
    public FixedSizeIterable<GenotypeRecord> getSampleGenotypeRecords(GeneticVariant variant) {
        final VcfRecord vcfRecord = getVcfRecord(variant);

        return new FixedSizeIterable<GenotypeRecord>() {
            @Override
            public Iterator<GenotypeRecord> iterator() {
                Iterable<VcfSample> samples = vcfRecord.getSamples();
                return Iterators.transform(samples.iterator(), new Function<VcfSample, GenotypeRecord>() {
                    @Override
                    public GenotypeRecord apply(VcfSample vcfSample) {
                        return toGenotypeRecord(vcfRecord, vcfSample);
                    }
                });
            }

            @Override
            public int size() {
                return vcfRecord.getNrSamples();
            }
        };
    }

    @Override
    public List<String> getSeqNames() {
        return new ArrayList<String>(tabixIndex.getSeqNames());
    }

    @Override
    public Iterable<GeneticVariant> getSequenceGeneticVariants(String seqName) {
        return getVariantsByRange(seqName, 0, Integer.MAX_VALUE);
    }

    @Override
    public Iterable<GeneticVariant> getVariantsByPos(final String seqName, final int startPos) {
        return getVariantsByRange(seqName, startPos - 1, startPos);
    }

    @Override
    public GeneticVariant getSnpVariantByPos(String seqName, int startPos) {

        //NOTE this special version will complete the iteration and therefor the file connection will be closed
        //This is kind of an ugly hack but is functional
        Iterable<GeneticVariant> variants = getVariantsByPos(seqName, startPos);
        GeneticVariant snp = null;
        for (GeneticVariant variant : variants) {
            if (snp == null && variant.isSnp()) {
                snp = variant;
            }
        }

        return snp;

    }

    @Override
    public Iterable<GeneticVariant> getVariantsByRange(final String seqName, final int rangeStart, final int rangeEnd) {

        if (rangeStart < 0) {
            throw new GenotypeDataException("Illegal start pos for VCF variant query: " + rangeStart);
        }

        return new Iterable<GeneticVariant>() {
            @Override
            public Iterator<GeneticVariant> iterator() {

                try {

                    ++totalRandomAccessRequest;
                    ++currentlyOpenFileHandlers;

                    return new Iterator<GeneticVariant>() {
                        private final BlockCompressedInputStream stream = new BlockCompressedInputStream(bzipVcfFile);
                        private final TabixIterator it = tabixIndex.queryTabixIndex(seqName, rangeStart, rangeEnd, stream);
                        private String line = readFirst(it);

                        private String readFirst(TabixIterator it) {
                            if (it == null) {
                                try {
                                    stream.close();
                                } catch (IOException e) {
                                    throw new GenotypeDataException(e);
                                }
                                --currentlyOpenFileHandlers;
                                ++closedFileHandlers;
                                return null;
                            } else {
                                try {
                                    String firstLine = it.next();
                                    if (firstLine == null) {
                                        --currentlyOpenFileHandlers;
                                        ++closedFileHandlers;
                                        it.close();
                                    }
                                    return firstLine;
                                } catch (IOException e) {
                                    throw new GenotypeDataException(e);
                                }
                            }
                        }

                        @Override
                        public boolean hasNext() {
                            return line != null;
                        }

                        @Override
                        public GeneticVariant next() {
                            VcfRecord vcfRecord = new VcfRecord(vcfMeta, StringUtils.split(line, '\t'));
                            try {
                                line = it.next();
                                if (line == null) {
                                    --currentlyOpenFileHandlers;
                                    ++closedFileHandlers;
                                    it.close();
                                }
                            } catch (IOException e) {
                                throw new GenotypeDataException(e);
                            }
                            return toGeneticVariant(vcfRecord);
                        }

                        @Override
                        public void remove() {
                            throw new UnsupportedOperationException();
                        }
                    };
                } catch (FileNotFoundException e) {
                    if (e.getMessage().endsWith("(Too many open files)")) {
                        throw new GenotypeDataException("VCF reader trying to open more file connections than allowed by operating system. Currently open connections: " + currentlyOpenFileHandlers + " total opened: " + totalRandomAccessRequest + " total closed: " + closedFileHandlers, e);
                    } else {
                        throw new GenotypeDataException(e);
                    }
                } catch (IOException e) {
                    throw new GenotypeDataException(e);
                }
            }
        };
    }

    private VcfRecord getVcfRecord(GeneticVariant variant) {
        if (!variant.equals(cachedGeneticVariant)) {
            TabixIterator it;
            String line;
            BlockCompressedInputStream stream = null;
            ++totalRandomAccessRequest;
            ++currentlyOpenFileHandlers;

            try {
                stream = new BlockCompressedInputStream(bzipVcfFile);
                it = tabixIndex.queryTabixIndex(variant.getSequenceName(), variant.getStartPos() - 1, variant.getStartPos(), stream);
                while ((line = it.next()) != null) {
                    VcfRecord vcfRecord = new VcfRecord(vcfMeta, StringUtils.split(line, '\t'));
                    if (variant.equals(toGeneticVariant(vcfRecord))) {
                        cachedVcfRecord = vcfRecord;
                        cachedGeneticVariant = variant;
                        break;
                    }
                }
                stream.close();
            } catch (FileNotFoundException e) {
                if (e.getMessage().endsWith("(Too many open files)")) {
                    throw new GenotypeDataException("VCF reader trying to open more file connections than allowed by operating system. Currently open connections: " + currentlyOpenFileHandlers + " total opened: " + totalRandomAccessRequest + " total closed: " + closedFileHandlers, e);
                } else {
                    throw new GenotypeDataException(e);
                }
            } catch (IOException e) {
                throw new GenotypeDataException(e);
            } finally {
                --currentlyOpenFileHandlers;
                ++closedFileHandlers;
                IOUtils.closeQuietly(stream);
            }
        }
        return cachedVcfRecord;
    }

    /**
     * Convert VcfRecord to GeneticVariant
     *
     * @param vcfRecord
     * @return
     */
    private GeneticVariant toGeneticVariant(VcfRecord vcfRecord) {
        List<String> identifiers = vcfRecord.getIdentifiers();
        int pos = vcfRecord.getPosition();
        String sequenceName = vcfRecord.getChromosome();
        Allele refAllele = vcfRecord.getReferenceAllele();
        List<Allele> altAlleles = vcfRecord.getAlternateAlleles();

        Map<String, Object> annotationMap = new HashMap<String, Object>();
        for (VcfInfo vcfInfo : vcfRecord.getInformation()) {
            annotationMap.put(vcfInfo.getKey(), vcfInfo.getVal());
        }
        annotationMap.put("VCF_Filter", vcfRecord.getFilterStatus());
        annotationMap.put("VCF_Qual", vcfRecord.getQuality());

        Alleles alleles;
        if (altAlleles == null || altAlleles.isEmpty()) {
            alleles = Alleles.createAlleles(refAllele);
        } else {
            ArrayList<Allele> allelesList = new ArrayList<Allele>(altAlleles.size() + 1);
            allelesList.add(refAllele);
            allelesList.addAll(altAlleles);
            alleles = Alleles.createAlleles(allelesList);
        }

        GeneticVariantMeta geneticVariantMeta = new VcfGeneticVariantMeta(vcfMeta, Arrays.asList(vcfRecord.getFormat()));
        GeneticVariant geneticVariant = ReadOnlyGeneticVariant.createVariant(geneticVariantMeta, identifiers, pos, sequenceName, annotationMap, variantProvider, alleles, refAllele);

        cachedVcfRecord = vcfRecord;
        cachedGeneticVariant = geneticVariant;
        return geneticVariant;
    }

    /**
     * Convert VcfSample to GenotypeRecord
     *
     * @param vcfRecord
     * @param vcfSample
     * @return
     */
    private GenotypeRecord toGenotypeRecord(VcfRecord vcfRecord, VcfSample vcfSample) {
        return new VcfGenotypeRecord(vcfMeta, vcfRecord, vcfSample);
    }

    /**
     * Getter for the preferred VCF genotype field to read from
     *
     * @return The preferred genotype field to read from.
     */
    public VcfGenotypeFormat getPreferredGenotypeField() {
        return genotypeFormatSupplier.getPreferredGenotypeFormat();
    }

    /**
     * Set the preferred VCF format to read from.
     * Fails if the requested vcf format has no included mapping in the this VcfGenotypeData
     *
     * @param preferredGenotypeFormat the preferred VCF format to read from.
     */
    public void setPreferredGenotypeFormat(String preferredGenotypeFormat) {
        this.setPreferredGenotypeFormat(
                new VcfGenotypeFormatSupplier(preferredGenotypeFormat != null ?
                VcfGenotypeFormat.valueOf(preferredGenotypeFormat) : null));
    }

    public void setPreferredGenotypeFormat(String preferredGenotypeFormat, boolean raiseExceptionIfUnavailable) {
        this.setPreferredGenotypeFormat(
                new VcfGenotypeFormatSupplier(preferredGenotypeFormat != null ?
                        VcfGenotypeFormat.valueOf(preferredGenotypeFormat) : null, raiseExceptionIfUnavailable));
    }

    public VcfGenotypeFormatSupplier getPreferredGenotypeFormat() {
        return genotypeFormatSupplier;
    }

    public void setPreferredGenotypeFormat(
            VcfGenotypeFormat preferredGenotypeField, String fieldIdentifier, boolean raiseExceptionIfUnavailable) {
        this.setPreferredGenotypeFormat(new VcfGenotypeFormatSupplier(
                preferredGenotypeField, fieldIdentifier,
                raiseExceptionIfUnavailable));
    }

    public void setPreferredGenotypeFormat(
            VcfGenotypeFormat preferredGenotypeField, String fieldIdentifier) {
        this.setPreferredGenotypeFormat(new VcfGenotypeFormatSupplier(
                preferredGenotypeField, fieldIdentifier));
    }

    public void setPreferredGenotypeFormat(VcfGenotypeFormatSupplier genotypeFormatSupplier) {
        this.genotypeFormatSupplier = genotypeFormatSupplier;
    }

}
