/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.downstreamer.gene;

import cern.colt.matrix.tdouble.DoubleMatrix2D;
import java.util.ArrayList;
import java.util.LinkedHashMap;

import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.LogManager;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.variant.GeneticVariant;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;

/**
 *
 * @author patri
 */
@Deprecated
public class GenotypeCorrelationGenotypes implements GenotypeCorrelationSource {

	private static final Logger LOGGER = LogManager.getLogger(GenotypeCorrelationGenotypes.class);
	private static final DoubleMatrixDataset<String, String> EMPTY_DATASET = new DoubleMatrixDataset<>(0, 0);

	private final RandomAccessGenotypeData referenceGenotypes;
	private final LinkedHashMap<String, Integer> sampleHash;

	public GenotypeCorrelationGenotypes(RandomAccessGenotypeData referenceGenotypes) {
		this.referenceGenotypes = referenceGenotypes;

		String[] samples = referenceGenotypes.getSampleNames();

		sampleHash = new LinkedHashMap<>(samples.length);

		int s = 0;
		for (String sample : samples) {
			sampleHash.put(sample, s++);
		}
	}

	@Override
	public DoubleMatrixDataset<String, String> getCorrelationMatrixForRange(String chr, int start, int stop) {

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
		GenePvalueCalculator.timeInLoadingGenotypeDosages += (timeStop - timeStart);

		LOGGER.debug(" * Variants found in locus: " + variantsDosages.size());

		if (dosageDataset.rows() == 0) {
			return EMPTY_DATASET;
		} else {
			timeStart = System.currentTimeMillis();
			DoubleMatrixDataset<String, String> corMatrix = dosageDataset.calculateCorrelationMatrix();
			timeStop = System.currentTimeMillis();
			GenePvalueCalculator.timeInCreatingGenotypeCorrelationMatrix += (timeStop - timeStart);
			return corMatrix;

		}

	}

}
