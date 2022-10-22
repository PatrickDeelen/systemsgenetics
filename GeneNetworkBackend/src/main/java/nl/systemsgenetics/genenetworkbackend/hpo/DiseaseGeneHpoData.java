/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.genenetworkbackend.hpo;

import cern.colt.matrix.tdouble.DoubleMatrix1D;
import com.opencsv.CSVParser;
import com.opencsv.CSVParserBuilder;
import com.opencsv.CSVReader;
import com.opencsv.CSVReaderBuilder;

import java.io.*;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Objects;
import java.util.Random;
import java.util.Set;
import java.util.function.Predicate;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;

import org.apache.commons.math3.stat.regression.SimpleRegression;
import org.apache.commons.math3.util.FastMath;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;

/**
 * @author patri
 */
public class DiseaseGeneHpoData {

	private final HashMap<String, HashSet<String>> geneToHpos;
	private final HashMap<String, HashSet<String>> diseaseToGenes;
	private final HashMap<DiseaseGene, HashSet<String>> diseaseGeneToHpos; // disease_gene
	private final SimpleRegression regression = new SimpleRegression();

	public DiseaseGeneHpoData(final File diseaseGeneHpoFile, HashMap<String, ArrayList<String>> ncbiToEnsgMap, HashMap<String, ArrayList<String>> hgncToEnsgMap, HashSet<String> exludedHpo, HashSet<String> includeGenes, String diseasePrefix) throws FileNotFoundException, IOException {

		geneToHpos = new HashMap<>();
		diseaseToGenes = new HashMap<>();
		diseaseGeneToHpos = new HashMap<>();

		Predicate<String> diseasePattern;
		if (diseasePrefix != null) {
			diseasePattern = Pattern.compile("^" + diseasePrefix).asPredicate();
		} else {
			diseasePattern = null;
		}

		final CSVParser hpoParser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();

		CSVReader hpoReader = null;
		if (diseaseGeneHpoFile.getName().endsWith(".gz")) {
			hpoReader = new CSVReaderBuilder(new BufferedReader(new InputStreamReader((new GZIPInputStream(new FileInputStream(diseaseGeneHpoFile)))))).withSkipLines(1).withCSVParser(hpoParser).build();
		} else {
			hpoReader = new CSVReaderBuilder(new BufferedReader(new FileReader(diseaseGeneHpoFile))).withSkipLines(1).withCSVParser(hpoParser).build();
		}

		String[] nextLine;
		while ((nextLine = hpoReader.readNext()) != null) {
			String disease = nextLine[0];
			String hgcnId = nextLine[1];
			String ncbiId = nextLine[2];
			String hpo = nextLine[3];

			if (diseasePattern != null && !diseasePattern.test(disease)) {
				continue;
			}

			if (exludedHpo != null && exludedHpo.contains(hpo)) {
				continue;
			}

			ArrayList<String> ensgIds = ncbiToEnsgMap.get(ncbiId);
			if (ensgIds == null) {
				ensgIds = hgncToEnsgMap.get(hgcnId);
			}
			if (ensgIds == null) {
				System.err.println("Missing mapping for gene: " + ncbiId + " " + hgcnId);
			} else if (ensgIds.size() > 1) {
				System.err.println("Skipping becasue multiple ENSG IDs for gene: " + ncbiId + " " + hgcnId);
			} else if (!includeGenes.contains(ensgIds.get(0))) {
				System.err.println("Skipping becasue gene not in include list: " + ncbiId + " " + hgcnId);
			} else {

				String ensgId = ensgIds.get(0);

				HashSet<String> geneHpos = geneToHpos.get(ensgId);
				if (geneHpos == null) {
					geneHpos = new HashSet<>();
					geneToHpos.put(ensgId, geneHpos);
				}

				geneHpos.add(hpo);

				HashSet<String> diseaseGenes = diseaseToGenes.get(disease);
				if (diseaseGenes == null) {
					diseaseGenes = new HashSet<>();
					diseaseToGenes.put(disease, diseaseGenes);
				}
				diseaseGenes.add(ensgId);

				DiseaseGene diseaseGene = new DiseaseGene(disease, ensgId);

				HashSet<String> diseaseGeneHpos = diseaseGeneToHpos.get(diseaseGene);
				if (diseaseGeneHpos == null) {
					diseaseGeneHpos = new HashSet<>();
					diseaseGeneToHpos.put(diseaseGene, diseaseGeneHpos);
				}
				diseaseGeneHpos.add(hpo);

			}

		}

	}

	public DiseaseGeneHpoData(HashMap<DiseaseGene, HashSet<String>> diseaseGeneToHpos) {

		this.diseaseGeneToHpos = diseaseGeneToHpos;

		geneToHpos = new HashMap<>();
		diseaseToGenes = new HashMap<>();

		for (Map.Entry<DiseaseGene, HashSet<String>> diseaseGeneToHposEntry : diseaseGeneToHpos.entrySet()) {

			DiseaseGene diseaseGene = diseaseGeneToHposEntry.getKey();
			HashSet<String> hpos = diseaseGeneToHposEntry.getValue();

			HashSet<String> geneHpos = geneToHpos.get(diseaseGene.getGene());
			if (geneHpos == null) {
				geneHpos = new HashSet<>();
				geneToHpos.put(diseaseGene.getGene(), geneHpos);
			}

			geneHpos.addAll(hpos);

			HashSet<String> diseaseGenes = diseaseToGenes.get(diseaseGene.getDisease());
			if (diseaseGenes == null) {
				diseaseGenes = new HashSet<>();
				diseaseToGenes.put(diseaseGene.getDisease(), diseaseGenes);
			}
			diseaseGenes.add(diseaseGene.getGene());

		}

	}

	/**
	 * Returns null if no phenotypes associated
	 *
	 * @param ensgId
	 * @return
	 */
	public Set<String> getEnsgHpos(String ensgId) {

		HashSet<String> geneHpos = geneToHpos.get(ensgId);

		if (geneHpos == null) {
			return null;
		} else {
			return Collections.unmodifiableSet(geneHpos);
		}

	}

	public Set<String> getDiseaseGenes() {
		return Collections.unmodifiableSet(geneToHpos.keySet());
	}

	public Set<String> getDiseases() {
		return Collections.unmodifiableSet(diseaseToGenes.keySet());
	}

	public Set<DiseaseGene> getDiseaseGeneHpos() {
		return Collections.unmodifiableSet(diseaseGeneToHpos.keySet());
	}

	/**
	 * Returns null if no disease genes are found
	 *
	 * @param disease
	 * @return
	 */
	public Set<String> getGenesForDisease(String disease) {
		HashSet<String> diseaseGenes = diseaseToGenes.get(disease);

		if (diseaseGenes == null) {
			return null;
		} else {
			return Collections.unmodifiableSet(diseaseGenes);
		}
	}

	/**
	 * Returns null if no phenotypes associated
	 *
	 * @param diseaseGene disease_gene
	 * @return
	 */
	public Set<String> getDiseaseEnsgHpos(DiseaseGene diseaseGene) {

		HashSet<String> hpos = diseaseGeneToHpos.get(diseaseGene);

		if (hpos == null) {
			return null;
		} else {
			return Collections.unmodifiableSet(hpos);
		}

	}

	public DiseaseGeneHpoData getPermutation() {
		return getPermutation(new Random(), null, null, 0, null, 0);
	}

	public DiseaseGeneHpoData getPermutation(long seed) {
		return getPermutation(new Random(seed), null, null, 0, null, 0);
	}

	public DiseaseGeneHpoData getPermutation(long seed, ArrayList<String> backgroundGenes) {
		return getPermutation(new Random(seed), backgroundGenes, null, 0, null, 0);
	}

	public DiseaseGeneHpoData getPermutation(ArrayList<String> backgroundGenes) {
		return getPermutation(new Random(), backgroundGenes, null, 0, null, 0);
	}

	public DiseaseGeneHpoData getPermutation(long seed, ArrayList<String> backgroundGenes, DoubleMatrixDataset<String, String> predictionMatrixSignificantCorrelationMatrix, double minCorrelationTomatch) {
		return getPermutation(new Random(seed), backgroundGenes, predictionMatrixSignificantCorrelationMatrix, minCorrelationTomatch, null, 0);
	}

	public DiseaseGeneHpoData getPermutation(long seed, ArrayList<String> backgroundGenes, DoubleMatrixDataset<String, String> predictionMatrixSignificantCorrelationMatrix, double minCorrelationTomatch, DoubleMatrixDataset<String, String> predictionMatrixSignificant, double minCorrelationToMatchGenes) {
		return getPermutation(new Random(seed), backgroundGenes, predictionMatrixSignificantCorrelationMatrix, minCorrelationTomatch, predictionMatrixSignificant, minCorrelationToMatchGenes);
	}

	private DiseaseGeneHpoData getPermutation(Random random, ArrayList<String> backgroundGenes, DoubleMatrixDataset<String, String> predictionMatrixSignificantCorrelationMatrix, double minCorrelationToMatchTerms, DoubleMatrixDataset<String, String> predictionMatrixSignificant, double minCorrelationToMatchGenes) {

		if (backgroundGenes == null) {
			backgroundGenes = new ArrayList(geneToHpos.keySet());
		}

		HashMap<DiseaseGene, HashSet<String>> randomDiseaseGeneToHpos = new HashMap<>();

		for (Map.Entry<DiseaseGene, HashSet<String>> diseaseGeneToHposEntry : this.diseaseGeneToHpos.entrySet()) {

			DiseaseGene diseaseGene = diseaseGeneToHposEntry.getKey();
			HashSet<String> hpos = diseaseGeneToHposEntry.getValue();

			String disease = diseaseGene.getDisease();
			String gene = diseaseGene.getGene();

			if (predictionMatrixSignificant != null && !predictionMatrixSignificant.containsRow(gene)) {
				continue;
			}

			HashSet<String> knownGenesForDisease = this.diseaseToGenes.get(disease);

			String randomReplacementGene;
			DiseaseGene randomDiseaseGene = null;
			boolean hpoOverlap;
			boolean hpoCorrelated;
			boolean genePredictionsCorrelated;

			int i = 0;
			boolean noRandomFound = false;
			int randomElement = -1;

			findRandomMatch:
			do {

				if (i++ >= 500000) {
					System.err.println("No random match found");
					noRandomFound = true;
					break;
				}

				if (backgroundGenes.isEmpty()) {
					System.err.println("No background genes left");
					noRandomFound = true;
					break;
				}

				hpoOverlap = false;
				hpoCorrelated = false;
				genePredictionsCorrelated = false;

				randomElement = random.nextInt(backgroundGenes.size());

				randomReplacementGene = backgroundGenes.get(randomElement);
				randomDiseaseGene = new DiseaseGene(disease, randomReplacementGene);
				HashSet<String> knownHposForRandomGene = this.geneToHpos.get(randomReplacementGene);

				if (knownHposForRandomGene != null) {
					for (String hpo : hpos) {
						if (knownHposForRandomGene.contains(hpo)) {
							hpoOverlap = true;
							continue findRandomMatch;
						}
					}
				}

				if (predictionMatrixSignificantCorrelationMatrix != null && knownHposForRandomGene != null) {
					//if already hpo overlap no need to do this

					hposLoop:
					for (String hpo : hpos) {

						if (predictionMatrixSignificantCorrelationMatrix.containsCol(hpo)) {

							for (String randomHpo : knownHposForRandomGene) {

								if (predictionMatrixSignificantCorrelationMatrix.containsCol(randomHpo) && predictionMatrixSignificantCorrelationMatrix.getElement(hpo, randomHpo) >= minCorrelationToMatchTerms) {
									hpoCorrelated = true;
									continue findRandomMatch;
								}

							}
						}
					}

				}

				if (predictionMatrixSignificant != null) {

					if (!predictionMatrixSignificant.containsRow(randomReplacementGene)) {
						genePredictionsCorrelated = true;//put to true to force selecting other gene
						continue findRandomMatch;
					}

					DoubleMatrix1D realGenePredictions = predictionMatrixSignificant.getRow(gene);
					DoubleMatrix1D randomGenePredictions = predictionMatrixSignificant.getRow(randomReplacementGene);

					for (int j = 0; j < realGenePredictions.size(); j++) {
						regression.addData(realGenePredictions.get(j), randomGenePredictions.get(j));
					}

					genePredictionsCorrelated = FastMath.abs(regression.getR()) > minCorrelationToMatchGenes;

				}

			} while (genePredictionsCorrelated | hpoCorrelated | hpoOverlap | knownGenesForDisease.contains(randomReplacementGene) | randomDiseaseGeneToHpos.containsKey(randomDiseaseGene));
			//geneToHpos.keySet().contains(randomReplacementGene) 

			if (!noRandomFound) {
				//backgroundGenes.remove(randomElement);
				randomDiseaseGeneToHpos.put(randomDiseaseGene, hpos);
			}

		}

		return new DiseaseGeneHpoData(randomDiseaseGeneToHpos);

	}

	public class DiseaseGene {

		private final String disease;
		private final String gene;

		public DiseaseGene(String disease, String gene) {
			this.disease = disease;
			this.gene = gene;
		}

		public String getDisease() {
			return disease;
		}

		public String getGene() {
			return gene;
		}

		@Override
		public int hashCode() {
			int hash = 3;
			hash = 97 * hash + Objects.hashCode(this.disease);
			hash = 97 * hash + Objects.hashCode(this.gene);
			return hash;
		}

		@Override
		public boolean equals(Object obj) {
			if (this == obj) {
				return true;
			}
			if (obj == null) {
				return false;
			}
			if (getClass() != obj.getClass()) {
				return false;
			}
			final DiseaseGene other = (DiseaseGene) obj;
			if (!Objects.equals(this.disease, other.disease)) {
				return false;
			}
			if (!Objects.equals(this.gene, other.gene)) {
				return false;
			}
			return true;
		}

		@Override
		public String toString() {
			return disease + "_" + gene;
		}

	}

}
