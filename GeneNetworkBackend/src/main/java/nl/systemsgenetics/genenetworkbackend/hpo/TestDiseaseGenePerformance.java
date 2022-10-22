/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.genenetworkbackend.hpo;

import umcg.genetica.io.hpo.HpoOntology;
import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import com.opencsv.CSVParser;
import com.opencsv.CSVParserBuilder;
import com.opencsv.CSVReader;
import com.opencsv.CSVReaderBuilder;
import com.opencsv.CSVWriter;

import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.Map;
import java.util.Random;
import java.util.Set;
import java.util.zip.GZIPInputStream;

import static nl.systemsgenetics.genenetworkbackend.ConvertHpoToMatrix.loadHgncToEnsgMap;
import static nl.systemsgenetics.genenetworkbackend.ConvertHpoToMatrix.loadNcbiToEnsgMap;
import static nl.systemsgenetics.genenetworkbackend.div.CalculateGenePredictability.loadSignificantTerms;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.moment.Variance;
import org.apache.commons.math3.stat.descriptive.rank.Median;
import org.apache.commons.math3.stat.ranking.NaNStrategy;
import org.apache.commons.math3.stat.ranking.NaturalRanking;
import org.apache.commons.math3.stat.ranking.TiesStrategy;
import org.biojava.nbio.ontology.Ontology;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;

/**
 *
 * @author patri
 */
public class TestDiseaseGenePerformance {

	/**
	 * @param args the command line arguments
	 * @throws java.lang.Exception
	 */
	public static void main(String[] args) throws Exception {

		final File diseaseGeneHpoFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\HPO\\135\\ALL_SOURCES_ALL_FREQUENCIES_diseases_to_genes_to_phenotypes.txt");
		final File ncbiToEnsgMapFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\ensgNcbiId.txt");
		final File hgncToEnsgMapFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\ensgHgnc.txt");
		final File ensgSymbolMappingFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\ensgHgnc.txt");
		final File predictionMatrixFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\Data31995Genes05-12-2017\\PCA_01_02_2018\\predictions2\\hpo_predictions.txt.gz");
		final File predictionMatrixCorrelationFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\Data31995Genes05-12-2017\\PCA_01_02_2018\\predictions\\hpo_predictions_pathwayCorrelation.txt");
		final File significantTermsFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\Data31995Genes05-12-2017\\PCA_01_02_2018\\predictions2\\hpo_predictions_bonSigTerms.txt");
		final double correctedPCutoff = 0.05;
		final File hpoOboFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\HPO\\135\\hp.obo");
		final File hpoPredictionInfoFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\Data31995Genes05-12-2017\\PCA_01_02_2018\\predictions2\\hpo_predictions_auc_bonferroni.txt");
		final File hposToExcludeFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\hpoToExclude.txt");
		final File skewnessFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\Data31995Genes05-12-2017\\PCA_01_02_2018\\predictions\\skewnessSummary.txt");
		final boolean randomize = false;
		final File annotationMatrixFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\Data31995Genes05-12-2017\\PCA_01_02_2018\\PathwayMatrix\\ALL_SOURCES_ALL_FREQUENCIES_phenotype_to_genes.txt_matrix.txt.gz");
		//final File backgroundForRandomize = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\Data31995Genes05-12-2017\\PCA_01_02_2018\\PathwayMatrix\\Ensembl2Reactome_All_Levels.txt_genesInPathways.txt");
		//final File backgroundForRandomize = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\expressedReactomeGenes.txt");
		final boolean randomizeCustomBackground = false;
		final int limitHpoTermUsage = 5; //0 is no limit
		
		Map<String, String> ensgSymbolMapping = loadEnsgToHgnc(ensgSymbolMappingFile);

		
		final StringBuilder outputFilePath = new StringBuilder("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\hpoDiseaseBenchmark");
		final ArrayList<String> backgroundGenes;
		
		
		if(limitHpoTermUsage > 0){
			outputFilePath.append("hpoLimit").append(limitHpoTermUsage);
		}
		
		if (randomize) {

			if (randomizeCustomBackground) {
				System.err.println("First need to fix so ranking list contains all genes in background list");
				return;
//				backgroundGenes = loadBackgroundGenes(backgroundForRandomize);
//				outputFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\hpoDiseaseBenchmarkRandomizedCustomBackground.txt");
			} else {
				backgroundGenes = null;
				outputFilePath.append("Randomized");
			}

		} else {
			backgroundGenes = null;
		}
		
		
		outputFilePath.append(".txt");
		final File outputFile = new File(outputFilePath.toString());

		final HashMap<String, ArrayList<String>> ncbiToEnsgMap = loadNcbiToEnsgMap(ncbiToEnsgMapFile);
		final HashMap<String, ArrayList<String>> hgncToEnsgMap = loadHgncToEnsgMap(hgncToEnsgMapFile);
		final HashSet<String> exludedHpo = loadHpoExclude(hposToExcludeFile);

		final SkewnessInfo skewnessInfo = new SkewnessInfo(skewnessFile);

		LinkedHashSet<String> significantTerms = loadSignificantTerms(significantTermsFile);

		DoubleMatrixDataset<String, String> predictionMatrix = DoubleMatrixDataset.loadDoubleData(predictionMatrixFile.getAbsolutePath());
		DoubleMatrixDataset<String, String> predictionMatrixSignificant = predictionMatrix.viewColSelection(significantTerms);

		DoubleMatrixDataset<String, String> predictionMatrixSignificantCorrelationMatrix = DoubleMatrixDataset.loadDoubleData(predictionMatrixCorrelationFile.getAbsolutePath());

		DiseaseGeneHpoData diseaseGeneHpoData = new DiseaseGeneHpoData(diseaseGeneHpoFile, ncbiToEnsgMap, hgncToEnsgMap, exludedHpo, new HashSet(predictionMatrix.getHashRows().keySet()), "OMIM");
		
		//NOTE if one would use a differnt background this needs to be updated
		HashSet<String> diseaseGenes = new HashSet<>(diseaseGeneHpoData.getDiseaseGenes());
		
		if (randomize) {
			diseaseGeneHpoData = diseaseGeneHpoData.getPermutation(1, backgroundGenes);
		}

		

		for (String gene : diseaseGenes) {
			if (!predictionMatrixSignificant.containsRow(gene)) {
				throw new Exception("Error: " + gene);
			}
		}

		int[] mapGeneIndexToDiseaseGeneIndex = new int[predictionMatrix.rows()];
		ArrayList<String> predictedGenes = predictionMatrix.getRowObjects();

		int g2 = 0;
		for (int g = 0; g < predictedGenes.size(); ++g) {
			mapGeneIndexToDiseaseGeneIndex[g] = diseaseGenes.contains(predictedGenes.get(g)) ? g2++ : -1;
		}

		DoubleMatrixDataset<String, String> annotationnMatrix = DoubleMatrixDataset.loadDoubleData(annotationMatrixFile.getAbsolutePath());
		DoubleMatrixDataset<String, String> annotationMatrixSignificant = annotationnMatrix.viewColSelection(significantTerms);

		HashMap<String, MeanSd> hpoMeanSds = calculatePathayMeansOfAnnotatedGenes(predictionMatrixSignificant, annotationMatrixSignificant);

		Map<String, PredictionInfo> predictionInfo = HpoFinder.loadPredictionInfo(hpoPredictionInfoFile);

		Ontology hpoOntology = HpoOntology.loadHpoOntology(hpoOboFile);

		HpoFinder hpoFinder = new HpoFinder(hpoOntology, predictionInfo);

		final int totalGenes = predictionMatrixSignificant.rows();
		final int totalDiseaseGenes = diseaseGenes.size();
		final double[] geneScores = new double[totalGenes];
		final double[] geneScoresDiseaseGenes = new double[totalDiseaseGenes];
		final NaturalRanking naturalRanking = new NaturalRanking(NaNStrategy.FAILED, TiesStrategy.MAXIMUM);

		CSVWriter writer = new CSVWriter(new FileWriter(outputFile), '\t', '\0', '\0', "\n");

		String[] outputLine = new String[16];
		int c = 0;
		outputLine[c++] = "Disease";
		outputLine[c++] = "Gene";
		outputLine[c++] = "Hgnc";
		outputLine[c++] = "Rank";
		outputLine[c++] = "RankAmongDiseaseGenes";
		outputLine[c++] = "Z-score";
		outputLine[c++] = "HPO_skewness";
		outputLine[c++] = "Other_mean_skewness";
		outputLine[c++] = "Other_max_skewness";
		outputLine[c++] = "HPO_phenotypic_match_score";
		outputLine[c++] = "HPO_count";
		outputLine[c++] = "HPO_sum_auc";
		outputLine[c++] = "HPO_mean_auc";
		outputLine[c++] = "HPO_median_auc";
		outputLine[c++] = "HPO_terms";
		outputLine[c++] = "HPO_terms_match_score";
		writer.writeNext(outputLine);

		Random random = new Random(1);

		Mean meanCalculator = new Mean();
		Median medianCalculator = new Median();

		for (DiseaseGeneHpoData.DiseaseGene diseaseGene : diseaseGeneHpoData.getDiseaseGeneHpos()) {

			String gene = diseaseGene.getGene();
			String disease = diseaseGene.getDisease();

			if (!predictionMatrixSignificant.containsRow(gene)) {
				continue;
			}

			Set<String> geneHpos = diseaseGeneHpoData.getDiseaseEnsgHpos(diseaseGene);

			LinkedHashSet<String> geneHposPredictable = new LinkedHashSet<>();

			for (String hpo : geneHpos) {
				geneHposPredictable.addAll(hpoFinder.getTermsToNames(hpoFinder.getPredictableTerms(hpo, correctedPCutoff)));
			}

			if (geneHposPredictable.isEmpty()) {
				continue;
			}

			if(limitHpoTermUsage > 0 && geneHposPredictable.size() > limitHpoTermUsage){
				ArrayList<String> hposTmp = new ArrayList(geneHposPredictable);
				geneHposPredictable.clear();
				
				for(int i = 0; i < limitHpoTermUsage ; ++i){
					geneHposPredictable.add(hposTmp.remove(random.nextInt(hposTmp.size())));
				}
				
			}
			
//			if(geneHposPredictable.size() > 1){
//				String hpoSelected = geneHposPredictable.toArray(new String[geneHposPredictable.size()])[random.nextInt(geneHposPredictable.size())];
//				geneHposPredictable = new LinkedHashSet<>(1);
//				geneHposPredictable.add(hpoSelected);
//			}
			DoubleMatrixDataset<String, String> predictionCaseTerms = predictionMatrixSignificant.viewColSelection(geneHposPredictable);
			DoubleMatrix2D predictionCaseTermsMatrix = predictionCaseTerms.getMatrix();

			double denominator = Math.sqrt(geneHposPredictable.size());

			for (int g = 0; g < totalGenes; ++g) {
				geneScores[g] = predictionCaseTermsMatrix.viewRow(g).zSum() / denominator;
				if (Double.isNaN(geneScores[g])) {
					geneScores[g] = 0;
				}

				g2 = mapGeneIndexToDiseaseGeneIndex[g];
				if (g2 >= 0) {
					geneScoresDiseaseGenes[g2] = geneScores[g];
				}

			}

			double[] geneRanks = naturalRanking.rank(geneScores);
			int diseaseGeneIndex = predictionMatrixSignificant.getRowIndex(gene);

			double[] geneRanksDiseaseGenes = naturalRanking.rank(geneScoresDiseaseGenes);
			int diseaseGeneIndexInDiseaseGenesOnly = mapGeneIndexToDiseaseGeneIndex[diseaseGeneIndex];

			double zscore = geneScores[diseaseGeneIndex];
			double rank = (totalGenes - geneRanks[diseaseGeneIndex]) + 1;
			double rankAmongDiseaseGenes = (totalDiseaseGenes - geneRanksDiseaseGenes[diseaseGeneIndexInDiseaseGenesOnly]) + 1;

			double hpoPhenotypicMatchScore = 0;
			StringBuilder individualMatchScore = new StringBuilder();
			boolean notFirst = false;
			int usedHpos = 0;
			
			double[] aucs = new double[geneHposPredictable.size()];
			double sumAucs = 0;
			
			int i = 0;
			for (String hpo : geneHposPredictable) {

				usedHpos++;

				MeanSd hpoMeanSd = hpoMeanSds.get(hpo);

				double hpoPredictionZ = predictionMatrixSignificant.getElement(gene, hpo);

				double hpoPredictionOutlierScore = ((hpoPredictionZ - hpoMeanSd.getMean()) / hpoMeanSd.getSd());

				if (notFirst) {
					individualMatchScore.append(';');
				}
				notFirst = true;

				individualMatchScore.append(hpoPredictionOutlierScore);

				hpoPhenotypicMatchScore += hpoPredictionOutlierScore;
				
				aucs[i++] = predictionInfo.get(hpo).getAuc();
				sumAucs += predictionInfo.get(hpo).getAuc();

			}
			
			double meanAuc = meanCalculator.evaluate(aucs);
			double medianAuc = medianCalculator.evaluate(aucs);

			if (usedHpos == 0) {
				hpoPhenotypicMatchScore = Double.NaN;
			} else {
				hpoPhenotypicMatchScore = hpoPhenotypicMatchScore / usedHpos;
			}

			String symbol = ensgSymbolMapping.get(gene);
			if (symbol == null) {
				symbol = "";
			}

			c = 0;
			outputLine[c++] = disease;
			outputLine[c++] = gene;
			outputLine[c++] = symbol;
			outputLine[c++] = String.valueOf(rank);
			outputLine[c++] = String.valueOf(rankAmongDiseaseGenes);
			outputLine[c++] = String.valueOf(zscore);
			outputLine[c++] = String.valueOf(skewnessInfo.getHpoSkewness(gene));
			outputLine[c++] = String.valueOf(skewnessInfo.getMeanSkewnessExHpo(gene));
			outputLine[c++] = String.valueOf(skewnessInfo.getMaxSkewnessExHpo(gene));
			outputLine[c++] = String.valueOf(hpoPhenotypicMatchScore);
			outputLine[c++] = String.valueOf(geneHposPredictable.size());
			outputLine[c++] = String.valueOf(sumAucs);
			outputLine[c++] = String.valueOf(meanAuc);
			outputLine[c++] = String.valueOf(medianAuc);
			outputLine[c++] = String.join(";", geneHposPredictable);
			outputLine[c++] = individualMatchScore.toString();
			writer.writeNext(outputLine);

		}

		writer.close();

	}

	public static HashSet<String> loadHpoExclude(File hposToExclude) throws IOException {

		LinkedHashSet<String> hpos = new LinkedHashSet<>();

		BufferedReader reader = new BufferedReader(new FileReader(hposToExclude));

		String line;
		while ((line = reader.readLine()) != null) {
			hpos.add(line);
		}

		return hpos;

	}

	private static HashMap<String, MeanSd> calculatePathayMeansOfAnnotatedGenes(DoubleMatrixDataset<String, String> predictionMatrixSignificant, DoubleMatrixDataset<String, String> annotationMatrixSignificant) {

		HashMap<String, MeanSd> pathwayMeanSdMap = new HashMap<>(predictionMatrixSignificant.columns());

		Mean meanCalculator = new Mean();
		Variance varianceCalculator = new Variance();

		for (String pathway : predictionMatrixSignificant.getColObjects()) {

			meanCalculator.clear();
			varianceCalculator.clear();

			DoubleMatrix1D pathwayPredictions = predictionMatrixSignificant.getCol(pathway);
			DoubleMatrix1D pathwayAnnotations = annotationMatrixSignificant.getCol(pathway);

			for (int g = 0; g < pathwayPredictions.size(); ++g) {
				if (pathwayAnnotations.get(g) != 0) {
					meanCalculator.increment(pathwayPredictions.getQuick(g));
					varianceCalculator.increment(pathwayPredictions.getQuick(g));
				}
			}

			double v = varianceCalculator.getResult();

			pathwayMeanSdMap.put(pathway, new MeanSd(meanCalculator.getResult(), v * v));

		}

		return pathwayMeanSdMap;

	}

	private static class MeanSd {

		private final double mean;
		private final double sd;

		public MeanSd(double mean, double sd) {
			this.mean = mean;
			this.sd = sd;
		}

		public double getMean() {
			return mean;
		}

		public double getSd() {
			return sd;
		}

	}

	public static ArrayList<String> loadBackgroundGenes(File backgroundForRandomize) throws IOException {

		ArrayList<String> genes = new ArrayList<>();

		BufferedReader reader = new BufferedReader(new FileReader(backgroundForRandomize));

		String line;
		while ((line = reader.readLine()) != null) {
			genes.add(line);
		}

		return genes;

	}

	private static Map<String, String> loadEnsgToHgnc(File mappingFile) throws IOException {

		final CSVParser parser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
		CSVReader reader = null;
		if (mappingFile.getName().endsWith(".gz")) {
			reader = new CSVReaderBuilder(new BufferedReader(new InputStreamReader((new GZIPInputStream(new FileInputStream(mappingFile)))))).withSkipLines(1).withCSVParser(parser).build();
		} else {
			reader = new CSVReaderBuilder(new BufferedReader(new FileReader(mappingFile))).withSkipLines(1).withCSVParser(parser).build();
		}

		HashMap<String, String> mapping = new HashMap<>();

		String[] nextLine;
		while ((nextLine = reader.readNext()) != null) {

			mapping.put(nextLine[0], nextLine[1]);

		}

		return Collections.unmodifiableMap(mapping);

	}

}
