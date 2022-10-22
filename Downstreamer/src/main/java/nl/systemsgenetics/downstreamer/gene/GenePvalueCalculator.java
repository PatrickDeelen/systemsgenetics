/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.downstreamer.gene;

import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tdouble.algo.decomposition.DenseDoubleEigenvalueDecomposition;
import cern.jet.math.tdouble.DoubleFunctions;
import ch.unil.genescore.vegas.Farebrother;
import com.opencsv.CSVParser;
import com.opencsv.CSVParserBuilder;
import com.opencsv.CSVReader;
import com.opencsv.CSVReaderBuilder;
import com.opencsv.CSVWriter;

import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.concurrent.ThreadLocalRandom;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.IntStream;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import me.tongfei.progressbar.ProgressBar;
import me.tongfei.progressbar.ProgressBarStyle;

import static nl.systemsgenetics.downstreamer.Downstreamer.formatMsForLog;

import nl.systemsgenetics.downstreamer.Downstreamer;
import nl.systemsgenetics.downstreamer.DownstreamerOptions;
import nl.systemsgenetics.downstreamer.io.IoUtils;
import org.apache.commons.math3.distribution.ChiSquaredDistribution;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.apache.log4j.Logger;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.variant.GeneticVariant;
import org.molgenis.genotype.variantFilter.VariantIdIncludeFilter;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import umcg.genetica.math.matrix2.DoubleMatrixDatasetFastSubsetLoader;
import umcg.genetica.math.matrix2.DoubleMatrixDatasetRowIterable;
import umcg.genetica.math.stats.PearsonRToZscoreBinned;
import umcg.genetica.math.stats.ZScores;

/**
 * @author patri
 */
public class GenePvalueCalculator {

	public static final int MAX_ROUND_1_RESCUE = 10000000;

	private static final Logger LOGGER = Logger.getLogger(GenePvalueCalculator.class);
	private static final DoubleMatrixDataset<String, String> EMPTY_DATASET = new DoubleMatrixDataset<>(0, 0);
	private static final int PERMUTATION_STEP = 1000;
	private static final int MIN_PERMUTATIONS = 10000;
	private static final double MIN_PVALUE_FAREBROTHER = 1e-12;
	private final int numberRandomPhenotypes;
	//private static final int NUMBER_PERMUTATION_NULL_GWAS = 10000;
	//private static final double NUMBER_PERMUTATION_NULL_GWAS_PLUS_1 = NUMBER_PERMUTATION_NULL_GWAS + 1;

	protected static long timeInCreatingGenotypeCorrelationMatrix = 0;
	protected static long timeInLoadingGenotypeDosages = 0;

	private static long totalTimeInThread = 0;
	private static long timeInPermutations = 0;
	private static long timeInPruningGenotypeCorrelationMatrix = 0;
	private static long timeInDoingPca = 0;
	private static long timeInLoadingZscoreMatrix = 0;
	private static long timeInCalculatingPvalue = 0;
	private static long timeInCalculatingRealSumChi2 = 0;
	private static long timeInComparingRealChi2ToPermutationChi2 = 0;
	private static long timeInNormalizingGenotypeDosages = 0;

	private static int countRanPermutationsForGene = 0;
	private static int countBasedPvalueOnPermutations = 0;
	private static int countBasedPvalueOnFarebrother = 0;
	private static int countBasedPvalueOnPermutationsAfterFailedFarebrother = 0;
	private static int countUseChi2DistForPvalue = 0;
	private static int countNoVariants = 0;

	private final RandomAccessGenotypeData referenceGenotypes;
	private final List<Gene> genes;
	private final int windowExtend;
	private final double maxR;
	private final int maxNrPermutations;
	private final int maxNrPermutationsRescue1;
	private final long maxNrPermutationsRescue2;
	//private final double nrPermutationsPlus1Double;
	private final String outputBasePath;
	private final double[] randomChi2;
	private final LinkedHashMap<String, Integer> sampleHash;//In order of the samples in the genotype data
	private final DoubleMatrixDataset<String, String> randomNormalizedPhenotypes;//Rows samples order of sample hash, cols phenotypes
	private final PearsonRToZscoreBinned r2zScore;//Based on number of samples with genotypes
	private final DoubleMatrixDataset<String, String> genePvalues;
	private final DoubleMatrixDataset<String, String> geneMaxSnpZscore;
	private final DoubleMatrixDataset<String, String> genePvaluesNullGwas;
	private final DoubleMatrixDataset<String, String> geneMaxSnpZscoreNullGwas;
	private final DoubleMatrixDataset<String, String> geneVariantCount;
	private final DoubleMatrixDataset<String, String> geneMaxPermutationCount;
	private final DoubleMatrixDataset<String, String> geneRuntime;
	private final int numberRealPheno;
	private final DoubleMatrixDatasetFastSubsetLoader geneVariantPhenotypeMatrixRowLoader;
	private final int[] genePValueDistributionPermuations;
	private final int[] genePValueDistributionChi2Dist;
	private final ProgressBar pb;
	//	private final double minPvaluePermutations;
//	private final double minPvaluePermutations2;
	private final boolean correctForLambdaInflation;
	private final double[] lambdaInflations;
	private final File debugFolder;
	private final HashMap<String, GeneticVariant> variantIdMap;
	private final HashMap<String, HashSet<String>> variantGeneMapping;
	private final CSVWriter variantPerGeneWriter;
	private final String[] variantPerGeneOutputLine;

	/**
	 * @param variantPhenotypeZscoreMatrixPath Binary matrix. rows: variants,
	 *                                         cols: phenotypes. Variants must have IDs identical to
	 *                                         genotypeCovarianceSource
	 * @param referenceGenotypes               Source data used to calculate correlations
	 *                                         between SNPs
	 * @param genes                            genes for which to calculate p-values.
	 * @param windowExtend                     number of bases to add left and right of gene window
	 * @param maxR                             max correlation between variants to use
	 * @param nrPermutations
	 * @param outputBasePath
	 * @param randomChi2
	 * @param correctForLambdaInflation
	 * @param nrSampleToUseForCorrelation
	 * @param nrSamplesToUseForNullBetas
	 * @throws java.io.IOException
	 */
	@SuppressWarnings("CallToThreadStartDuringObjectConstruction")
	public GenePvalueCalculator(String variantPhenotypeZscoreMatrixPath, RandomAccessGenotypeData referenceGenotypes, List<Gene> genes,
								int windowExtend, double maxR, int nrPermutations, long nrRescuePermutation, String outputBasePath,
								double[] randomChi2, boolean correctForLambdaInflation, final int nrSampleToUseForCorrelation,
								final int nrSamplesToUseForNullBetas, final File debugFolder,
								final File variantGeneMappingFile, final File usedVariantsPerGeneFile) throws IOException, Exception {

		this.referenceGenotypes = referenceGenotypes;
		this.genes = genes;
		this.windowExtend = windowExtend;
		this.maxR = maxR;
		this.maxNrPermutations = nrPermutations;
		this.maxNrPermutationsRescue1 = nrRescuePermutation < MAX_ROUND_1_RESCUE ? (int) nrRescuePermutation : MAX_ROUND_1_RESCUE;
		this.maxNrPermutationsRescue2 = nrRescuePermutation;
		this.debugFolder = debugFolder;

//		this.minPvaluePermutations = 0.5 / (this.maxNrPermutations2 + 1);
//		this.minPvaluePermutations2 = 0.5 / (this.maxNrPermutations2 + 1);
		//this.nrPermutationsPlus1Double = nrPermutations + 1;
		this.outputBasePath = outputBasePath;
		this.randomChi2 = randomChi2;
		this.correctForLambdaInflation = correctForLambdaInflation;
		this.numberRandomPhenotypes = nrSampleToUseForCorrelation + nrSamplesToUseForNullBetas;

		if (usedVariantsPerGeneFile != null) {
			if (usedVariantsPerGeneFile.getName().endsWith(".gz")) {
				variantPerGeneWriter = new CSVWriter(new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(usedVariantsPerGeneFile)))), '\t', '\0', '\0', "\n");
			} else {
				variantPerGeneWriter = new CSVWriter(new FileWriter(usedVariantsPerGeneFile), '\t', '\0', '\0', "\n");
			}
			variantPerGeneOutputLine = new String[2];
			int c = 0;
			variantPerGeneOutputLine[c++] = "Variant";
			variantPerGeneOutputLine[c++] = "Gene";
			variantPerGeneWriter.writeNext(variantPerGeneOutputLine);
		} else {
			variantPerGeneWriter = null;
			variantPerGeneOutputLine = null;
		}

		if (variantGeneMappingFile != null) {
			this.variantGeneMapping = loadVariantGeneMapping(variantGeneMappingFile);

			HashSet<String> variantInGeneMapping = new HashSet<>();
			for (HashSet<String> variants : variantGeneMapping.values()) {
				variantInGeneMapping.addAll(variants);
			}
			this.variantIdMap = this.referenceGenotypes.getVariantIdMap(new VariantIdIncludeFilter(variantInGeneMapping));

		} else {
			this.variantGeneMapping = null;
			this.variantIdMap = null;
		}

		if (nrPermutations % PERMUTATION_STEP != 0) {
			throw new Exception("Number of permutaitons must be dividable by: " + PERMUTATION_STEP);
		}

		String[] samples = referenceGenotypes.getSampleNames();
		sampleHash = new LinkedHashMap<>(samples.length);

		int s = 0;
		for (String sample : samples) {
			sampleHash.put(sample, s++);
		}

		randomNormalizedPhenotypes = generateRandomNormalizedPheno(sampleHash, numberRandomPhenotypes);
		r2zScore = new PearsonRToZscoreBinned(10000000, sampleHash.size());//10000000
		String variantPhenotypeZscoreMatrixPathColFile = variantPhenotypeZscoreMatrixPath + ".cols.txt";
		if (!new File(variantPhenotypeZscoreMatrixPathColFile).canRead()) {
			variantPhenotypeZscoreMatrixPathColFile += ".gz";
			if (!new File(variantPhenotypeZscoreMatrixPathColFile).canRead()) {
				throw new FileNotFoundException("Could not find file: " + variantPhenotypeZscoreMatrixPath + ".cols.txt or " + variantPhenotypeZscoreMatrixPath + ".cols.txt.gz");
			}
		}
		final List<String> phenotypes = IoUtils.readMatrixAnnotations(new File(variantPhenotypeZscoreMatrixPathColFile));

		numberRealPheno = phenotypes.size();
		final int numberGenes = genes.size();

		if (this.correctForLambdaInflation) {

			this.lambdaInflations = new double[numberRealPheno];
			final double y = new ChiSquaredDistribution(1).inverseCumulativeProbability(0.5);
			DoubleMatrixDatasetRowIterable gwasStreamer = new DoubleMatrixDatasetRowIterable(variantPhenotypeZscoreMatrixPath);
			DescriptiveStatistics[] medianCalculators = new DescriptiveStatistics[numberRealPheno];
			for (int p = 0; p < numberRealPheno; ++p) {
				medianCalculators[p] = new DescriptiveStatistics(gwasStreamer.getNrRows());
			}

			for (double[] variantZscores : gwasStreamer) {
				for (int p = 0; p < numberRealPheno; ++p) {
					medianCalculators[p].addValue(variantZscores[p] * variantZscores[p]);
				}
			}

			IntStream.range(0, numberRealPheno).parallel().forEach(p -> {
				final double medianChi2 = medianCalculators[p].getPercentile(50);
				lambdaInflations[p] = medianChi2 / y;
				//if(LOGGER.isDebugEnabled()){
				LOGGER.info("Pheno_" + phenotypes.get(p) + " median chi2: " + medianChi2 + " lambda inflation: " + lambdaInflations[p]);
				//}
			});
		} else {
			this.lambdaInflations = null;
		}

		// Result matrix. Rows: genes, Cols: phenotypes
		genePvalues = new DoubleMatrixDataset<>(createGeneHashRows(genes), createHashColsFromList(phenotypes));
		geneMaxSnpZscore = new DoubleMatrixDataset<>(createGeneHashRows(genes), createHashColsFromList(phenotypes));

		//Result matrix null GWAS. Rows: genes, Cols: phenotypes
		genePvaluesNullGwas = new DoubleMatrixDataset<>(createGeneHashRows(genes), randomNormalizedPhenotypes.getHashCols());
		geneMaxSnpZscoreNullGwas = new DoubleMatrixDataset<>(createGeneHashRows(genes), randomNormalizedPhenotypes.getHashCols());

		geneVariantCount = new DoubleMatrixDataset<>(createGeneHashRows(genes), createHashColsFromList(Arrays.asList(new String[]{"count"})));
		geneMaxPermutationCount = new DoubleMatrixDataset<>(createGeneHashRows(genes), createHashColsFromList(Arrays.asList(new String[]{"maxPermutations"})));
		geneRuntime = new DoubleMatrixDataset<>(createGeneHashRows(genes), createHashColsFromList(Arrays.asList(new String[]{"runtime"})));

		geneVariantPhenotypeMatrixRowLoader = new DoubleMatrixDatasetFastSubsetLoader(variantPhenotypeZscoreMatrixPath);

		genePValueDistributionPermuations = new int[21];//used to create histogram 
		genePValueDistributionChi2Dist = new int[21];//used to create histogram 

		pb = new ProgressBar("Gene p-value calculations", numberGenes, ProgressBarStyle.ASCII);

		//All genes are indipendant
		// IntStream.range(0, numberGenes).parallel().forEach((int geneI) -> {
//			for (int geneI = 0; geneI < numberGenes; ++geneI) {
		long startThread = System.currentTimeMillis();

		//Results are writen in genePvalues
		//runGene(geneI);
		final AtomicInteger count = new AtomicInteger(0);
		List<Thread> threads = new ArrayList<>(DownstreamerOptions.getNumberOfThreadsToUse());
		Downstreamer.ThreadErrorHandler threadErrorHandler = new Downstreamer.ThreadErrorHandler("gene p-value calculator");

		for (int i = 0; i < DownstreamerOptions.getNumberOfThreadsToUse(); ++i) {

			Thread worker = new Thread(new CalculatorThread(count));
			worker.start();
			worker.setUncaughtExceptionHandler(threadErrorHandler);
			threads.add(worker);

		}

		boolean running;
		do {
			running = false;
			for (Thread thread : threads) {
				if (thread.isAlive()) {
					running = true;
				}
			}

			try {
				Thread.sleep(500);
			} catch (InterruptedException ex) {
			}

		} while (running);

		long endThread = System.currentTimeMillis();

		totalTimeInThread += (endThread - startThread);

		pb.close();
		if (variantPerGeneWriter != null) {
			variantPerGeneWriter.close();
		}

		LOGGER.info("countRanPermutationsForGene: " + countRanPermutationsForGene);
		LOGGER.info("countBasedPvalueOnPermutations: " + countBasedPvalueOnPermutations);
		LOGGER.info("countBasedPvalueOnPermutationsAfterFailedFarebrother: " + countBasedPvalueOnPermutationsAfterFailedFarebrother);
		LOGGER.info("countBasedPvalueOnFarebrother: " + countBasedPvalueOnFarebrother);
		LOGGER.info("countUseChi2DistForPvalue: " + countUseChi2DistForPvalue);
		LOGGER.info("countNoVariants: " + countNoVariants);

		LOGGER.info("timeInLoadingGenotypeDosages: " + formatMsForLog(timeInLoadingGenotypeDosages));
		LOGGER.info("timeInNormalizingGenotypeDosages: " + formatMsForLog(timeInNormalizingGenotypeDosages));
		LOGGER.info("timeInCreatingGenotypeCorrelationMatrix: " + formatMsForLog(timeInCreatingGenotypeCorrelationMatrix));
		LOGGER.info("timeInPruningGenotypeCorrelationMatrix: " + formatMsForLog(timeInPruningGenotypeCorrelationMatrix));
		LOGGER.info("timeInDoingPca: " + formatMsForLog(timeInDoingPca));
		LOGGER.info("timeInPermutation: " + formatMsForLog(timeInPermutations));
		LOGGER.info("timeInLoadingZscoreMatrix: " + formatMsForLog(timeInLoadingZscoreMatrix));
		LOGGER.info("timeInCalculatingRealSumChi2: " + formatMsForLog(timeInCalculatingRealSumChi2));
		LOGGER.info("timeInComparingRealChi2ToPermutationChi2: " + formatMsForLog(timeInComparingRealChi2ToPermutationChi2));
		LOGGER.info("timeInCalculatingPvalue: " + formatMsForLog(timeInCalculatingPvalue));
		LOGGER.info("totalTimeInThread: " + formatMsForLog(totalTimeInThread));

		LOGGER.info("-----------------------");
		LOGGER.info("Gene p-value histrogram chi2 dist");
		for (double histCount : genePValueDistributionChi2Dist) {
			LOGGER.info(histCount);
		}

		LOGGER.info("Gene p-value histrogram permuations");
		for (double histCount : genePValueDistributionPermuations) {
			LOGGER.info(histCount);
		}

		LOGGER.info("-----------------------");

	}

	/**
	 * @return gene p-value matrix for each phenotype. rows: genes in same order
	 * as genes list, cols: phenotypes
	 */
	public DoubleMatrixDataset<String, String> getGenePvalues() {
		return genePvalues;
	}

	public DoubleMatrixDataset<String, String> getGeneMaxSnpZscore() {
		return geneMaxSnpZscore;
	}

	public DoubleMatrixDataset<String, String> getGenePvaluesNullGwas() {
		return genePvaluesNullGwas;
	}

	public DoubleMatrixDataset<String, String> getGeneMaxSnpZscoreNullGwas() {
		return geneMaxSnpZscoreNullGwas;
	}

	public DoubleMatrixDataset<String, String> getGeneVariantCount() {
		return geneVariantCount;
	}

	public DoubleMatrixDataset<String, String> getGeneMaxPermutationCount() {
		return geneMaxPermutationCount;
	}

	public DoubleMatrixDataset<String, String> getGeneRuntime() {
		return geneRuntime;
	}

	private void runGene(int geneI, final double[] geneChi2SumNull) throws IOException, Exception {

		long geneTimeStart = System.currentTimeMillis();

		long timeStart;
		long timeStop;

		Gene gene = genes.get(geneI);

		//Rows: samples, cols: variants
		final DoubleMatrixDataset<String, String> variantScaledDosages;
		final DoubleMatrixDataset<String, String> variantScaledDosagesPruned;
		final DoubleMatrixDataset<String, String> variantCorrelationsPruned;
		final int variantCorrelationsPrunedRows;

		{

			final DoubleMatrixDataset<String, String> variantScaledDosagesRange;
			if (windowExtend >= 0) {
				DoubleMatrixDataset<String, String> tmp = loadVariantScaledDosageMatrix(gene.getChr(), gene.getStart() - windowExtend, gene.getEnd() + windowExtend);

				if (tmp.columns() == 0) {
					//no variants in range
					variantScaledDosagesRange = null;
				} else {
					variantScaledDosagesRange = tmp;
				}

			} else {
				variantScaledDosagesRange = null;
			}

			final DoubleMatrixDataset<String, String> variantScaledDosagesIds;
			if (variantGeneMapping != null && variantGeneMapping.containsKey(gene.getGene())) {
				DoubleMatrixDataset<String, String> tmp = loadVariantScaledDosageMatrix(variantGeneMapping.get(gene.getGene()));
				if (tmp.columns() == 0) {
					//no variants in range
					variantScaledDosagesIds = null;
				} else {
					variantScaledDosagesIds = tmp;
				}
			} else {
				variantScaledDosagesIds = null;
			}

			if (variantScaledDosagesRange == null && variantScaledDosagesIds == null) {
				variantScaledDosages = new DoubleMatrixDataset<>(0, 0);
			} else if (variantScaledDosagesRange == null) {
				variantScaledDosages = variantScaledDosagesIds;
			} else if (variantScaledDosagesIds == null) {
				variantScaledDosages = variantScaledDosagesRange;
			} else {
				//Both are not null so they need to be merged

				//First identify variants in de ID matrix that are not in the range matrix
				ArrayList<String> variantsFoundBasedOnId = variantScaledDosagesIds.getColObjects();
				Iterator<String> variantsFoundBasedOnIdIt = variantsFoundBasedOnId.iterator();
				while (variantsFoundBasedOnIdIt.hasNext()) {
					String variantId = variantsFoundBasedOnIdIt.next();
					if (variantScaledDosagesRange.containsCol(variantId)) {
						variantsFoundBasedOnIdIt.remove();
					}
				}

				//variantsFoundBasedOnId now only contains variants not yet in variantScaledDosagesRange
				if (variantsFoundBasedOnId.isEmpty()) {
					variantScaledDosages = variantScaledDosagesRange;
				} else {
					variantScaledDosages = DoubleMatrixDataset.concatColumns(variantScaledDosagesRange, variantScaledDosagesIds.viewColSelection(variantsFoundBasedOnId));
				}

			}

			// Calculate the correlation matrix between variants (LD matrix sigma)
			timeStart = System.currentTimeMillis();
			final DoubleMatrixDataset<String, String> variantCorrelations = variantScaledDosages.calculateCorrelationMatrixOnNormalizedColumns();
			timeStop = System.currentTimeMillis();
			timeInCreatingGenotypeCorrelationMatrix += (timeStop - timeStart);

			// Save correlation matrix
			if (LOGGER.isDebugEnabled() & variantCorrelations.rows() > 1) {
				variantCorrelations.save(new File(debugFolder, gene.getGene() + "_variantCorMatrix.txt.gz"));
			}

			// Write which variants are linked to which gene
			if (variantPerGeneWriter != null) {
				int c;
				synchronized (variantPerGeneWriter) {
					for (String variant : variantScaledDosages.getHashCols().keySet()) {
						c = 0;
						variantPerGeneOutputLine[c++] = variant;
						variantPerGeneOutputLine[c++] = gene.getGene();
						variantPerGeneWriter.writeNext(variantPerGeneOutputLine);
					}
					variantPerGeneWriter.flush();
				}
			}

			// Prune the correlation matrix (sigma) for genes with strong correlation as described by PASCAL
			timeStart = System.currentTimeMillis();
			variantCorrelationsPruned = pruneCorrelationMatrix(variantCorrelations, maxR);
			variantCorrelationsPrunedRows = variantCorrelationsPruned.rows();
			timeStop = System.currentTimeMillis();
			timeInPruningGenotypeCorrelationMatrix += (timeStop - timeStart);

			// Final pruned correlation matrix
			variantScaledDosagesPruned = variantScaledDosages.viewColSelection(variantCorrelationsPruned.getHashCols().keySet());
		}

		// Define null GWAS z-score by correlating a gausian dists (randomNormalizedPhenotypes)to the SNPs
		// in that region (variantScaledDosagesPruned)
		final DoubleMatrixDataset<String, String> nullGwasZscores;
		if (variantCorrelationsPrunedRows > 0) {
			// nullGwasZscores will first contain peason r values but this will be converted to Z-score using a lookup table;
			nullGwasZscores = DoubleMatrixDataset.correlateColumnsOf2ColumnNormalizedDatasets(randomNormalizedPhenotypes, variantScaledDosagesPruned);
			r2zScore.inplaceRToZ(nullGwasZscores);
		} else {
			nullGwasZscores = EMPTY_DATASET;
		}

		// Save the pruned correlation matrix to disk
		if (LOGGER.isDebugEnabled() & variantCorrelationsPrunedRows > 1) {
			variantCorrelationsPruned.save(new File(debugFolder, gene.getGene() + "_variantCorMatrixPruned.txt.gz"));
		}
		timeStop = System.currentTimeMillis();
		timeInCreatingGenotypeCorrelationMatrix += (timeStop - timeStart);

		// Count the number of pruned variants per gene
		geneVariantCount.setElementQuick(geneI, 0, variantCorrelationsPrunedRows);

		int currentNumberPermutationsCalculated;
		final ThreadLocalRandom rnd = ThreadLocalRandom.current();
		final double[] lambdas;
		final int lambdasLength;
		if (variantCorrelationsPrunedRows > 1) {

			// Determine the eigen decomposition of the SNP correlation matrix
			timeStart = System.currentTimeMillis();
			final DenseDoubleEigenvalueDecomposition eig = new DenseDoubleEigenvalueDecomposition(variantCorrelationsPruned.getMatrix());
			//final Jama.EigenvalueDecomposition eig = eigenValueDecomposition(variantCorrelationsPruned.getMatrixAs2dDoubleArray());

			// Get the eigenvalues for each eigenvector
			final DoubleMatrix1D eigenValues = eig.getRealEigenvalues().viewFlip();
			final long eigenValuesLenght = eigenValues.size();

			// Method below is from PASCAL to select relevant eigen values
			// <NOTE OBB> I didn't see anything in the paper about this, altough I kinda get why (I think)
			// <NOTE OBB> You don't want to include zero eigenvalues here. But on the other hand, does that matter
			// <NOTE OBB> as the eigenvalue approaches zero lambda * chi-sqr will just be zero, and add very little to
			// <NOTE OBB> the sum.
			double sumPosEigen = 0;
			for (int i = 0; i < eigenValuesLenght; i++) {
				double e = eigenValues.getQuick(i);
				if (e > 0) {
					sumPosEigen += e;
				}
			}

			// <NOTE OBB> why is it harcoded at 10k??
			final double cutoff = sumPosEigen / 10000; // Only use components that explain significant part of variation

			int eigenValuesToUse = 0;
			for (int i = 0; i < eigenValuesLenght; i++) {
				sumPosEigen -= eigenValues.getQuick(i);
				eigenValuesToUse++;
				if (sumPosEigen < cutoff) {
					break;
				}
			}

			lambdas = eigenValues.viewPart(0, eigenValuesToUse).toArray();
			lambdasLength = eigenValuesToUse;

			// Save the significant eigenvalues
			if (LOGGER.isDebugEnabled()) {
				saveEigenValues(lambdas, new File(debugFolder, gene.getGene() + "_eigenValues.txt"));
			}

			// Logging of running times
			timeStop = System.currentTimeMillis();
			timeInDoingPca += (timeStop - timeStart);
			timeStart = System.currentTimeMillis();
			timeStop = System.currentTimeMillis();
			timeInPermutations += (timeStop - timeStart);
			countRanPermutationsForGene++;

			// Generate a distribution of gaussian chi-sqr values weighted by the eigenvalues (lambdas)
			// Bit unclear coding, results is saved in geneChi2SumNull as this array is re-cycled.
			runPermutationsUsingEigenValues(geneChi2SumNull, lambdas, randomChi2, 0, MIN_PERMUTATIONS, rnd, lambdasLength);
			currentNumberPermutationsCalculated = MIN_PERMUTATIONS;

		} else {
			lambdas = null;
			lambdasLength = 0;
			currentNumberPermutationsCalculated = 0;
		}

		timeStart = System.currentTimeMillis();

		// Load current variants from variantPhenotypeMatrix
		final DoubleMatrixDataset<String, String> geneVariantPhenotypeMatrix;
		synchronized (this) {
			geneVariantPhenotypeMatrix = geneVariantPhenotypeMatrixRowLoader.loadSubsetOfRowsBinaryDoubleData(variantCorrelationsPruned.getHashRows().keySet());
		}

		if (LOGGER.isDebugEnabled()) {
			geneVariantPhenotypeMatrix.save(new File(debugFolder, gene.getGene() + "_variantPvalues.txt.gz"));
		}

		// <NOTE OBB> Why not just correct for lambda inflation once?
		if (correctForLambdaInflation) {
			final DoubleMatrix2D geneVariantPhenotypeMatrixInternal = geneVariantPhenotypeMatrix.getMatrix();
			for (int v = 0; v < geneVariantPhenotypeMatrix.rows(); ++v) {
				for (int p = 0; p < geneVariantPhenotypeMatrix.columns(); ++p) {
					geneVariantPhenotypeMatrixInternal.setQuick(v, p, Math.sqrt((geneVariantPhenotypeMatrixInternal.getQuick(v, p) * geneVariantPhenotypeMatrixInternal.getQuick(v, p)) / lambdaInflations[p]));
				}
			}
			if (LOGGER.isDebugEnabled()) {
				geneVariantPhenotypeMatrix.save(new File(debugFolder, gene.getGene() + "_variantPvaluesLambdaCorrected.txt.gz"));
			}
		}

		timeStop = System.currentTimeMillis();
		timeInLoadingZscoreMatrix += (timeStop - timeStart);

		Farebrother farebrother = null;

		for (int phenoI = 0; phenoI < numberRealPheno; ++phenoI) {

			if (variantCorrelationsPrunedRows > 1) {

				//PASCAL would now check if there is only a single significant lambda. For us this cannot hapend because of the pruning of SNPs
				timeStart = System.currentTimeMillis();

				DoubleMatrix1D phenoPvalues = geneVariantPhenotypeMatrix.getCol(phenoI);

				final double geneChi2Sum = phenoPvalues.aggregate(DoubleFunctions.plus, DoubleFunctions.square);
				final double geneMaxZscore = phenoPvalues.aggregate(DoubleFunctions.max, DoubleFunctions.abs);

				geneMaxSnpZscore.setElementQuick(geneI, phenoI, geneMaxZscore);

				timeStop = System.currentTimeMillis();
				timeInCalculatingRealSumChi2 += (timeStop - timeStart);

				timeStart = System.currentTimeMillis();

				int countNullLargerChi2ThanReal = 0;
				int perm = 0;
				long currentNumberPermutationsForThisPhenoGeneCombo = MIN_PERMUTATIONS - PERMUTATION_STEP;

				do {

					//LOGGER.debug("Start do with current permutations " + currentNumberPermutations + " x = " + x + " end: " + (currentNumberPermutations + PERMUTATION_STEP) + "?" + ((currentNumberPermutations + PERMUTATION_STEP) < maxNrPermutations));
					if (currentNumberPermutationsForThisPhenoGeneCombo >= currentNumberPermutationsCalculated) {
						runPermutationsUsingEigenValues(geneChi2SumNull, lambdas, randomChi2, currentNumberPermutationsCalculated, currentNumberPermutationsCalculated + PERMUTATION_STEP, rnd, lambdasLength);
						currentNumberPermutationsCalculated += PERMUTATION_STEP;
					}
					currentNumberPermutationsForThisPhenoGeneCombo += PERMUTATION_STEP;

					while (perm < currentNumberPermutationsForThisPhenoGeneCombo) {
						if (geneChi2Sum < geneChi2SumNull[perm++]) {
							countNullLargerChi2ThanReal++;
						}
					}

				} while (countNullLargerChi2ThanReal < 20 && currentNumberPermutationsForThisPhenoGeneCombo < maxNrPermutations);

//				if (currentNumberPermutationsForThisTraitGene == maxNrPermutations2) {
//					File nullFile = new File(outputBasePath + "_" + gene.getGene() + "_nullDistribution_" + variantCorrelationsPrunedRows + "variants" + "_xIs" + x + ".txt");
//
//					BufferedWriter nullWriter = new BufferedWriter(new FileWriter(nullFile));
//					for (int i = 0; i < geneChi2SumNull.length; ++i) {
//						nullWriter.write(String.valueOf(geneChi2SumNull[i]));
//						nullWriter.write('\n');
//					}
//					nullWriter.close();
//				}
				timeStop = System.currentTimeMillis();
				timeInComparingRealChi2ToPermutationChi2 += (timeStop - timeStart);

				timeStart = System.currentTimeMillis();

				final double p;

				LOGGER.debug("Done initial permutation for:\t" + gene.getGene() + "\tphenoI:\t" + phenoI + "\tcountNullLargerChi2ThanReal:\t" + countNullLargerChi2ThanReal + "\tused permutations:\t" + currentNumberPermutationsForThisPhenoGeneCombo);

				if (countNullLargerChi2ThanReal < 20) {

					//permutations not able to estimate p-value
					if (farebrother == null) {
						farebrother = new Farebrother(lambdas);
					}
					double farebrotherP = farebrother.probQsupx(geneChi2Sum);

					if (farebrother.getIfault() != 0) {
						//farebrother failed best we can do is some more permutation p-value

						//First fast permutations that can be used for other phenotypes / permutations
						do {
							//LOGGER.debug("Start do with current permutations " + currentNumberPermutations + " x = " + x + " end: " + (currentNumberPermutations + PERMUTATION_STEP) + "?" + ((currentNumberPermutations + PERMUTATION_STEP) < maxNrPermutations));
							if (currentNumberPermutationsForThisPhenoGeneCombo >= currentNumberPermutationsCalculated) {
								runPermutationsUsingEigenValues(geneChi2SumNull, lambdas, randomChi2, currentNumberPermutationsCalculated, currentNumberPermutationsCalculated + PERMUTATION_STEP, rnd, lambdasLength);
								currentNumberPermutationsCalculated += PERMUTATION_STEP;
							}
							currentNumberPermutationsForThisPhenoGeneCombo += PERMUTATION_STEP;

							while (perm < currentNumberPermutationsForThisPhenoGeneCombo) {
								if (geneChi2Sum < geneChi2SumNull[perm++]) {
									countNullLargerChi2ThanReal++;
								}
							}

						} while (countNullLargerChi2ThanReal < 20 && currentNumberPermutationsForThisPhenoGeneCombo < maxNrPermutationsRescue1);

						LOGGER.debug("Done first rescue permutation for:\t" + gene.getGene() + "\tphenoI:\t" + phenoI + "\tcountNullLargerChi2ThanReal:\t" + countNullLargerChi2ThanReal + "\tused permutations:\t" + currentNumberPermutationsForThisPhenoGeneCombo);

						//Fall back to slower purmtations that are not saved
						while (countNullLargerChi2ThanReal < 5 && currentNumberPermutationsForThisPhenoGeneCombo < maxNrPermutationsRescue2) {
							double weightedChi2Perm = 0;
							for (int g = -1; ++g < lambdasLength; ) {
								//double randomZ = rnd.nextGaussian();

								//Hopefully more efficient nextGaussian because no need to save second value from pair as double object
								double v1, v2, s;
								do {
									v1 = 2 * rnd.nextDouble() - 1; // between -1 and 1
									v2 = 2 * rnd.nextDouble() - 1; // between -1 and 1
									s = v1 * v1 + v2 * v2;
								} while (s >= 1 || s == 0);
								double multiplier = StrictMath.sqrt(-2 * StrictMath.log(s) / s);
								double randomZ = v1 * multiplier;
								weightedChi2Perm += lambdas[g] * (randomZ * randomZ);

								if (++g < lambdasLength) {
									randomZ = v2 * multiplier;
									weightedChi2Perm += lambdas[g] * (randomZ * randomZ);
								} else {
									break;//Will happen automaticly but would require second check
								}

							}
							if (geneChi2Sum < weightedChi2Perm) {
								countNullLargerChi2ThanReal++;
							}
							currentNumberPermutationsForThisPhenoGeneCombo++;
						}
						p = (countNullLargerChi2ThanReal + 0.5) / (double) (currentNumberPermutationsForThisPhenoGeneCombo + 1);
						countBasedPvalueOnPermutationsAfterFailedFarebrother++;
						LOGGER.debug(gene.getGene() + " permutation: " + currentNumberPermutationsForThisPhenoGeneCombo + " count null larger than real: " + countNullLargerChi2ThanReal + " pvalue:" + p);

					} else {
						//We noticed farebrother p-values below 1e-12 followed a very strange distribution
						p = farebrotherP <= 1e-12 ? 1e-12 : farebrotherP;
						LOGGER.debug("Using farebrother pvalue for:\t" + gene.getGene() + "\tphenoI:\t" + phenoI + "\tfarebrotherP:\t" + p);
						countBasedPvalueOnFarebrother++;
					}
				} else {
					countBasedPvalueOnPermutations++;
					p = (countNullLargerChi2ThanReal + 0.5) / (double) (currentNumberPermutationsForThisPhenoGeneCombo + 1);
					LOGGER.debug(gene.getGene() + " permutation: " + currentNumberPermutationsForThisPhenoGeneCombo + " count null larger than real: " + countNullLargerChi2ThanReal + " pvalue:" + p);
				}

				genePValueDistributionPermuations[(int) (20d * p)]++;
				genePvalues.setElementQuick(geneI, phenoI, p);

				timeStop = System.currentTimeMillis();
				timeInCalculatingPvalue += (timeStop - timeStart);

			} else if (variantCorrelationsPrunedRows == 1) {

				//Always row 0
				double p = ZScores.zToP(-Math.abs(geneVariantPhenotypeMatrix.getElementQuick(0, phenoI)));
				if (p == 1) {
					p = 0.99999d;
				}
				if (p < MIN_PVALUE_FAREBROTHER) {
					p = MIN_PVALUE_FAREBROTHER;
				}

				genePValueDistributionChi2Dist[(int) (20d * p)]++;
				genePvalues.setElementQuick(geneI, phenoI, p);

				countUseChi2DistForPvalue++;

			} else {
				//no variants in or near gene
				//genePValueDistribution[(int) (20d * 0.99999d)]++;
				genePvalues.setElementQuick(geneI, phenoI, Double.NaN);
				countNoVariants++;
			}

		}

		geneMaxPermutationCount.setElementQuick(geneI, 0, currentNumberPermutationsCalculated);

		int nullUsingFarebrother = 0;
		int nullUsingRescue1 = 0;
		int nullUsingRescue2 = 0;

		// Do exactly the same thing but now for the null GWAS
		for (int nullPhenoI = 0; nullPhenoI < numberRandomPhenotypes; ++nullPhenoI) {

			if (variantCorrelationsPrunedRows > 1) {

				timeStart = System.currentTimeMillis();

				DoubleMatrix1D phenoPvalues = nullGwasZscores.getCol(nullPhenoI);

				final double geneChi2Sum = phenoPvalues.aggregate(DoubleFunctions.plus, DoubleFunctions.square);
				final double geneMaxZscore = phenoPvalues.aggregate(DoubleFunctions.max, DoubleFunctions.abs);

				geneMaxSnpZscoreNullGwas.setElementQuick(geneI, nullPhenoI, geneMaxZscore);

				timeStop = System.currentTimeMillis();
				timeInCalculatingRealSumChi2 += (timeStop - timeStart);

				timeStart = System.currentTimeMillis();

				int countNullLargerChi2ThanReal = 0;
				int perm = 0;
				long currentNumberPermutationsForThisPhenoGeneCombo = MIN_PERMUTATIONS - PERMUTATION_STEP;

				do {

					//LOGGER.debug("Start do with current permutations " + currentNumberPermutations + " x = " + x + " end: " + (currentNumberPermutations + PERMUTATION_STEP) + "?" + ((currentNumberPermutations + PERMUTATION_STEP) < maxNrPermutations));
					if (currentNumberPermutationsForThisPhenoGeneCombo >= currentNumberPermutationsCalculated) {
						runPermutationsUsingEigenValues(geneChi2SumNull, lambdas, randomChi2, currentNumberPermutationsCalculated, currentNumberPermutationsCalculated + PERMUTATION_STEP, rnd, lambdasLength);
						currentNumberPermutationsCalculated += PERMUTATION_STEP;
					}
					currentNumberPermutationsForThisPhenoGeneCombo += PERMUTATION_STEP;

					while (perm < currentNumberPermutationsForThisPhenoGeneCombo) {
						if (geneChi2Sum < geneChi2SumNull[perm++]) {
							countNullLargerChi2ThanReal++;
						}
					}

				} while (countNullLargerChi2ThanReal < 20 && currentNumberPermutationsForThisPhenoGeneCombo < maxNrPermutations);

				timeStop = System.currentTimeMillis();
				timeInComparingRealChi2ToPermutationChi2 += (timeStop - timeStart);

				timeStart = System.currentTimeMillis();

				final double p;

				//For permutation phenotypes use threshold of 5 not 20
				if (countNullLargerChi2ThanReal < 5) {
					//permutations not able to estimate p-value
					if (farebrother == null) {
						farebrother = new Farebrother(lambdas);
					}

					++nullUsingFarebrother;

					double farebrotherP = farebrother.probQsupx(geneChi2Sum);

					if (farebrother.getIfault() != 0) {
						//farebrother failed best we can do is some more permutation p-value

						nullUsingRescue1++;

						//First fast permutations that can be used for other phenotypes / permutations
						do {
							//LOGGER.debug("Start do with current permutations " + currentNumberPermutations + " x = " + x + " end: " + (currentNumberPermutations + PERMUTATION_STEP) + "?" + ((currentNumberPermutations + PERMUTATION_STEP) < maxNrPermutations));
							if (currentNumberPermutationsForThisPhenoGeneCombo >= currentNumberPermutationsCalculated) {
								runPermutationsUsingEigenValues(geneChi2SumNull, lambdas, randomChi2, currentNumberPermutationsCalculated, currentNumberPermutationsCalculated + PERMUTATION_STEP, rnd, lambdasLength);
								currentNumberPermutationsCalculated += PERMUTATION_STEP;
							}
							currentNumberPermutationsForThisPhenoGeneCombo += PERMUTATION_STEP;

							while (perm < currentNumberPermutationsForThisPhenoGeneCombo) {
								if (geneChi2Sum < geneChi2SumNull[perm++]) {
									countNullLargerChi2ThanReal++;
								}
							}

						} while (countNullLargerChi2ThanReal < 5 && currentNumberPermutationsForThisPhenoGeneCombo < maxNrPermutationsRescue1);

						if (countNullLargerChi2ThanReal < 5) {
							++nullUsingRescue2;
						}

						//Fall back to slower purmtations that are not saved
						while (countNullLargerChi2ThanReal < 5 && currentNumberPermutationsForThisPhenoGeneCombo < maxNrPermutationsRescue2) {
							double weightedChi2Perm = 0;
							for (int g = -1; ++g < lambdasLength; ) {
								//double randomZ = rnd.nextGaussian();

								//Hopefully more efficient nextGaussian because no need to save second value from pair as double object
								double v1, v2, s;
								do {
									v1 = 2 * rnd.nextDouble() - 1; // between -1 and 1
									v2 = 2 * rnd.nextDouble() - 1; // between -1 and 1
									s = v1 * v1 + v2 * v2;
								} while (s >= 1 || s == 0);
								double multiplier = StrictMath.sqrt(-2 * StrictMath.log(s) / s);
								double randomZ = v1 * multiplier;
								weightedChi2Perm += lambdas[g] * (randomZ * randomZ);

								if (++g < lambdasLength) {
									randomZ = v2 * multiplier;
									weightedChi2Perm += lambdas[g] * (randomZ * randomZ);
								} else {
									break;//Will happen automaticly but would require second check
								}

							}
							if (geneChi2Sum < weightedChi2Perm) {
								countNullLargerChi2ThanReal++;
							}
							currentNumberPermutationsForThisPhenoGeneCombo++;
						}
						p = (countNullLargerChi2ThanReal + 0.5) / (double) (currentNumberPermutationsForThisPhenoGeneCombo + 1);
						//countBasedPvalueOnPermutationsAfterFailedFarebrother++;

					} else {
						//We noticed farebrother p-values below 1e-12 followed a very strange distribution
						p = farebrotherP <= 1e-12 ? 1e-12 : farebrotherP;
						//countBasedPvalueOnFarebrother++;
					}
				} else {
					//countBasedPvalueOnPermutations++;
					p = (countNullLargerChi2ThanReal + 0.5) / (double) (currentNumberPermutationsForThisPhenoGeneCombo + 1);

				}

				genePvaluesNullGwas.setElementQuick(geneI, nullPhenoI, p);

				timeStop = System.currentTimeMillis();
				timeInCalculatingPvalue += (timeStop - timeStart);

			} else if (variantCorrelationsPrunedRows == 1) {

				//Always row 0
				double p = ZScores.zToP(-Math.abs(nullGwasZscores.getElementQuick(0, nullPhenoI)));
				if (p == 1) {
					p = 0.99999d;
				}
				if (p < MIN_PVALUE_FAREBROTHER) {
					p = MIN_PVALUE_FAREBROTHER;
				}

				genePvaluesNullGwas.setElementQuick(geneI, nullPhenoI, p);

			} else {
				//no variants in or near gene
				//genePValueDistribution[(int) (20d * 0.99999d)]++;
				genePvaluesNullGwas.setElementQuick(geneI, nullPhenoI, Double.NaN);
			}

		}

		LOGGER.debug("Random phenotype permutations\t" + gene.getGene() + "\tnullUsingFarebrother\t" + nullUsingFarebrother + "\tnullUsingRescue1\t" + nullUsingRescue1 + "\tnullUsingRescue2\t" + nullUsingRescue2);

		geneRuntime.setElementQuick(geneI, 0, (System.currentTimeMillis() - geneTimeStart));

	}

	/**
	 * The dosages for each variants will be scaled to have mean of 0 and sd of
	 * 1. This will allow fast correlation calculations
	 *
	 * @param chr
	 * @param start
	 * @param stop
	 * @return
	 */
	private DoubleMatrixDataset<String, String> loadVariantScaledDosageMatrix(String chr, int start, int stop) {

		start = start < 0 ? 0 : start;

		LOGGER.debug("Query genotype data: " + chr + ":" + start + "-" + stop);

		//ArrayList<GeneticVariant> variants = new ArrayList<>(64);
		ArrayList<float[]> variantsDosages = new ArrayList<>(64);
		LinkedHashMap<String, Integer> variantHash = new LinkedHashMap<>(64);

		long timeStart = System.currentTimeMillis();

		int v = 0;
		synchronized (this) {
			for (GeneticVariant variant : referenceGenotypes.getVariantsByRange(chr, start, stop)) {
				//variants.add(variant);
				variantsDosages.add(variant.getSampleDosages());
				variantHash.put(variant.getPrimaryVariantId(), v++);
			}
		}

		DoubleMatrixDataset<String, String> dosageDataset = new DoubleMatrixDataset<>(sampleHash, variantHash);

		DoubleMatrix2D dosageMatrix = dosageDataset.getMatrix();

		v = 0;
		for (float[] variantDosages : variantsDosages) {
			for (int s = 0; s < variantDosages.length; ++s) {
				dosageMatrix.setQuick(s, v, variantDosages[s]);
			}
			v++;
		}

		long timeStop = System.currentTimeMillis();
		timeInLoadingGenotypeDosages += (timeStop - timeStart);

		LOGGER.debug(" * Variants found based on IDs: " + variantsDosages.size());

		timeStart = System.currentTimeMillis();

		//Inplace normalize per variants to mean of 0 and sd of 1 too 
		dosageDataset.normalizeColumns();

		timeStop = System.currentTimeMillis();
		timeInNormalizingGenotypeDosages += (timeStop - timeStart);

		return dosageDataset;

	}

	/**
	 * The dosages for each variants will be scaled to have mean of 0 and sd of
	 * 1. This will allow fast correlation calculations
	 *
	 * @param chr
	 * @param start
	 * @param stop
	 * @return
	 */
	private DoubleMatrixDataset<String, String> loadVariantScaledDosageMatrix(HashSet<String> variantsToLoad) {

		LOGGER.debug("Loading " + variantsToLoad.size() + "  variants based on IDs");

		//ArrayList<GeneticVariant> variants = new ArrayList<>(64);
		ArrayList<float[]> variantsDosages = new ArrayList<>(64);
		LinkedHashMap<String, Integer> variantHash = new LinkedHashMap<>(64);

		long timeStart = System.currentTimeMillis();

		int v = 0;
		synchronized (this) {
			for (String variantToLoad : variantsToLoad) {
				//variants.add(variant);
				GeneticVariant variant = variantIdMap.get(variantToLoad);
				if (variant != null) {
					variantsDosages.add(variant.getSampleDosages());
					variantHash.put(variant.getPrimaryVariantId(), v++);
				}
			}
		}

		DoubleMatrixDataset<String, String> dosageDataset = new DoubleMatrixDataset<>(sampleHash, variantHash);

		DoubleMatrix2D dosageMatrix = dosageDataset.getMatrix();

		v = 0;
		for (float[] variantDosages : variantsDosages) {
			for (int s = 0; s < variantDosages.length; ++s) {
				dosageMatrix.setQuick(s, v, variantDosages[s]);
			}
			v++;
		}

		long timeStop = System.currentTimeMillis();
		timeInLoadingGenotypeDosages += (timeStop - timeStart);

		LOGGER.debug(" * Variants found in locus: " + variantsDosages.size());

		timeStart = System.currentTimeMillis();

		//Inplace normalize per variants to mean of 0 and sd of 1 too 
		dosageDataset.normalizeColumns();

		timeStop = System.currentTimeMillis();
		timeInNormalizingGenotypeDosages += (timeStop - timeStart);

		return dosageDataset;

	}

	public static LinkedHashMap<String, Integer> createGeneHashRows(final List<Gene> genes) {
		LinkedHashMap<String, Integer> geneHashRows = new LinkedHashMap<>(genes.size());
		for (int geneI = 0; geneI < genes.size(); ++geneI) {
			geneHashRows.put(genes.get(geneI).getGene(), geneI);
		}
		return geneHashRows;
	}

	public static LinkedHashMap<String, Integer> createHashColsFromList(final List<String> list) {
		LinkedHashMap<String, Integer> hashCols = new LinkedHashMap<>(list.size());
		for (int i = 0; i < list.size(); ++i) {
			hashCols.put(list.get(i), i);
		}
		return hashCols;
	}

	/**
	 * Generate a distribution of sum chi-sqr values of nperm size, weighted by the eigenvalue lambda.
	 * For each permutation calculate sum(lambda * chi-sqr) over all values of lambda.
	 * <p>
	 * Fills a subset of the null
	 *
	 * @param geneChi2SumNull vector to fill will null chi2 values
	 * @param eigenValues
	 * @param randomChi2      reference distribution
	 * @param start           zero based count
	 * @param stop            zero based count, so must be < nrPermutation @return @param
	 *                        rnd @pa r am eigenValuesLengt h
	 */
	public static void runPermutationsUsingEigenValues(final double[] geneChi2SumNull, final double[] eigenValues, double[] randomChi2, final int start, final int stop, final ThreadLocalRandom rnd, final int eigenValuesLength) {

		final int randomChi2Size = randomChi2.length;

		// LOGGER.debug("start: " + start + " stop: " + stop);
		int i = start;
		for (int perm = start; perm < stop; perm++) {
			double weightedChi2Perm = 0;
			for (int g = 0; g < eigenValuesLength; g++) {

				if (i < randomChi2Size) {
					weightedChi2Perm += eigenValues[g] * randomChi2[i++];
				} else {
					double randomZ = rnd.nextGaussian();
					weightedChi2Perm += eigenValues[g] * (randomZ * randomZ);
				}

			}
			geneChi2SumNull[perm] = weightedChi2Perm;
		}

	}

	private static void saveEigenValues(double[] eigenValues, File file) throws IOException {
		final CSVWriter eigenWriter;
		if(file.getName().endsWith(".gz")){
			eigenWriter = new CSVWriter(new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(file)))), '\t', '\0', '\0', "\n");
		} else {
			eigenWriter = new CSVWriter(new FileWriter(file), '\t', '\0', '\0', "\n");
		}

		final String[] outputLine = new String[2];
		int c = 0;
		outputLine[c++] = "Component";
		outputLine[c++] = "EigenValue";
		eigenWriter.writeNext(outputLine);

		for (int i = 0; i < eigenValues.length; ++i) {

			c = 0;
			outputLine[c++] = "PC" + (i + 1);
			outputLine[c++] = String.valueOf(eigenValues[i]);
			eigenWriter.writeNext(outputLine);

		}

		eigenWriter.close();

	}

	protected static DoubleMatrixDataset<String, String> pruneCorrelationMatrix(DoubleMatrixDataset<String, String> correlationMatrixForRange, double maxR) {

		ArrayList<String> variantNames = correlationMatrixForRange.getRowObjects();
		LinkedHashSet<String> includedVariants = new LinkedHashSet<>(correlationMatrixForRange.rows());

		rows:
		for (int r = 0; r < correlationMatrixForRange.rows(); ++r) {
			cols:
			for (int c = 0; c < r; ++c) {
				if (Math.abs(correlationMatrixForRange.getElementQuick(r, c)) >= maxR && includedVariants.contains(variantNames.get(c))) {
					continue rows;
				}
			}
			includedVariants.add(variantNames.get(r));
		}

		LOGGER.debug(" * Variants after pruning high r: " + includedVariants.size());

		return correlationMatrixForRange.viewSelection(includedVariants, includedVariants);

	}

	private static DoubleMatrixDataset<String, String> generateRandomNormalizedPheno(LinkedHashMap<String, Integer> sampleHash, final int numberRandomPhenotypes) {

		LinkedHashMap<String, Integer> phenoHash = new LinkedHashMap<>();

		for (int i = 0; i < numberRandomPhenotypes; ++i) {
			phenoHash.put("RanPheno" + (i + 1), i);
		}

		DoubleMatrixDataset<String, String> phenoData = new DoubleMatrixDataset<>(sampleHash, phenoHash);

		DoubleMatrix2D phenoMatrix = phenoData.getMatrix();

		final int sampleCount = sampleHash.size();

		IntStream.range(0, numberRandomPhenotypes).parallel().forEach(pi -> {

			final ThreadLocalRandom rnd = ThreadLocalRandom.current();
			for (int s = 0; s < sampleCount; ++s) {
				phenoMatrix.setQuick(s, pi, rnd.nextGaussian());
			}

		});

		phenoData.normalizeColumns();

		return phenoData;

	}

	private static HashMap<String, HashSet<String>> loadVariantGeneMapping(File variantGeneMappingFile) throws FileNotFoundException, IOException {

		final CSVParser parser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();

		CSVReader reader = null;
		if (variantGeneMappingFile.getName().endsWith(".gz")) {
			reader = new CSVReaderBuilder((new BufferedReader(new InputStreamReader((new GZIPInputStream(new FileInputStream(variantGeneMappingFile))))))).withCSVParser(parser).build();
		} else {
			reader = new CSVReaderBuilder(new BufferedReader(new FileReader(variantGeneMappingFile))).withCSVParser(parser).build();
		}


		HashMap<String, HashSet<String>> geneVariantMapping = new HashMap<>();

		String[] nextLine;
		while ((nextLine = reader.readNext()) != null) {
			HashSet<String> geneVariant = geneVariantMapping.get(nextLine[1]);
			if (geneVariant == null) {
				geneVariant = new HashSet<>();
				geneVariantMapping.put(nextLine[1], geneVariant);
			}
			geneVariant.add(nextLine[0]);
		}

		return geneVariantMapping;

	}

	private class CalculatorThread implements Runnable {

		private final AtomicInteger counter;

		public CalculatorThread(AtomicInteger counter) {
			this.counter = counter;
		}

		@Override
		public void run() {

			int i;
			while ((i = counter.getAndIncrement()) < genes.size()) {
				try {
					final double[] geneChi2SumNull = new double[maxNrPermutationsRescue1];//The array will be recyceld. But the content will be overwritten
					runGene(i, geneChi2SumNull);
					pb.step();
				} catch (Exception ex) {
					throw new RuntimeException(ex);
				}
			}

		}

	}

}
