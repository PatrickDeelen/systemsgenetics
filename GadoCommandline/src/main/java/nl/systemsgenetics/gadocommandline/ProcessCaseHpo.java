/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.gadocommandline;

import com.opencsv.CSVParser;
import com.opencsv.CSVParserBuilder;
import com.opencsv.CSVReader;
import com.opencsv.CSVReaderBuilder;
import com.opencsv.CSVWriter;

import java.io.*;
import java.text.ParseException;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.zip.GZIPInputStream;

import org.apache.log4j.Logger;
import org.biojava.nbio.ontology.Ontology;
import org.biojava.nbio.ontology.Term;

/**
 * @author patri
 */
public class ProcessCaseHpo {

	private static final Logger LOGGER = Logger.getLogger(ProcessCaseHpo.class);

	/**
	 * @param options
	 * @throws IOException
	 * @throws FileNotFoundException
	 * @throws ParseException
	 */
	public static void process(GadoOptions options) throws IOException, FileNotFoundException, ParseException {

		final File hpoOboFile = options.getHpoOboFile();
		final File hpoPredictionInfoFile = options.getPredictionInfoFile();

		final File caseHpo = options.getInputCaseHpoFile();
		final File outputFile = new File(options.getOutputBasePath());

		final double correctedPCutoff = 0.05;

		Map<String, PredictionInfo> predictionInfo = HpoFinder.loadPredictionInfo(hpoPredictionInfoFile);

		Ontology hpoOntology = HpoFinder.loadHpoOntology(hpoOboFile);

		HpoFinder hpoFinder = new HpoFinder(hpoOntology, predictionInfo);

		//Map<String, String> updatedHpoId = loadUpdatedIds(updatedIdFile);

		final CSVParser parser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
		CSVReader reader = null;
		if (caseHpo.getName().endsWith(".gz")) {
			reader = new CSVReaderBuilder(new BufferedReader(
					new InputStreamReader((new GZIPInputStream(new FileInputStream(caseHpo))))
			)).withCSVParser(parser).build();
		} else {
			reader = new CSVReaderBuilder(new BufferedReader(new FileReader(caseHpo))).withCSVParser(parser).build();
		}


		CSVWriter writer = new CSVWriter(new FileWriter(outputFile), '\t', '\0', '\0', "\n");

		String[] outputLine = new String[6];
		int c = 0;
		outputLine[c++] = "Sample";
		outputLine[c++] = "SelectedHpo";
		outputLine[c++] = "SelectedHpoDescription";
		outputLine[c++] = "OriginalHpo";
		outputLine[c++] = "OriginalHpoDescription";
		outputLine[c++] = "ExcludeFromPrioritisation";
		writer.writeNext(outputLine);

		int sampleCounter = 0;

		String[] nextLine;
		while ((nextLine = reader.readNext()) != null) {

			++sampleCounter;

			String sampleId = nextLine[0];
			HashSet<String> sampleHpo = new HashSet<>();

			for (int i = 1; i < nextLine.length; i++) {

				String hpo = nextLine[i];

				if (hpo.length() == 0) {
					continue;
				}

//				if (updatedHpoId.containsKey(hpo)) {
//					hpo = updatedHpoId.get(hpo);
//				}

				if (sampleHpo.add(hpo)) {

					if (!hpoOntology.containsTerm(hpo)) {
						LOGGER.info("Warning term not found for: " + hpo + " of sample: " + sampleId);
					}

					Term hpoTerm = hpoOntology.getTerm(hpo);
					PredictionInfo info = predictionInfo.get(hpo);

					if (info == null || info.getCorrectedP() > correctedPCutoff) {
						//in case of no prediction or bad prediction

						List<Term> alternativeTerms = hpoFinder.getPredictableTerms(hpoTerm, correctedPCutoff);

						if (alternativeTerms.isEmpty()) {
							LOGGER.info("No alternative found for: " + hpo + " of sample: " + sampleId + " this term will be ignored");
						}

						for (Term alternativeTerm : alternativeTerms) {
							c = 0;
							outputLine[c++] = sampleId;
							outputLine[c++] = alternativeTerm.getName();
							outputLine[c++] = alternativeTerm.getDescription();
							outputLine[c++] = hpo;
							outputLine[c++] = hpoTerm.getDescription();
							outputLine[c++] = "";
							writer.writeNext(outputLine);
						}

					} else {

						c = 0;
						outputLine[c++] = sampleId;
						outputLine[c++] = hpo;
						outputLine[c++] = hpoTerm.getDescription();
						outputLine[c++] = "";
						outputLine[c++] = "";
						outputLine[c++] = "";
						writer.writeNext(outputLine);

					}

				}

			}

		}

		writer.close();

		LOGGER.info("Processed HPO terms of " + sampleCounter + " samples");

	}

	private static Map<String, String> loadUpdatedIds(File updatedIdFile) throws IOException {

		final CSVParser parser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
		CSVReader reader = null;
		if (updatedIdFile.getName().endsWith(".gz")) {
			reader = new CSVReaderBuilder(new BufferedReader(
					new InputStreamReader((new GZIPInputStream(new FileInputStream(updatedIdFile))))
			)).withCSVParser(parser).build();
		} else {
			reader = new CSVReaderBuilder(new BufferedReader(new FileReader(updatedIdFile))).withCSVParser(parser).build();
		}

		HashMap<String, String> updates = new HashMap<>();

		String[] nextLine;
		while ((nextLine = reader.readNext()) != null) {

			updates.put(nextLine[0], nextLine[1]);

		}

		return Collections.unmodifiableMap(updates);

	}

}
