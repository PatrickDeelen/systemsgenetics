package nl.systemsgenetics.genenetworkbackend;

import com.opencsv.CSVParser;
import com.opencsv.CSVParserBuilder;
import com.opencsv.CSVReader;
import com.opencsv.CSVReaderBuilder;

import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.Map;
import java.util.zip.GZIPInputStream;

import umcg.genetica.math.matrix2.DoubleMatrixDataset;

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 * @author patri
 */
public class ConvertReactomeToMatrix {

	/**
	 * @param args the command line arguments
	 * @throws java.io.IOException
	 * @throws java.lang.Exception
	 */
	public static void main(String[] args) throws IOException, Exception {

		final File pathwayFile = new File("D:\\UMCG\\Genetica\\Projects\\GeneNetwork\\Reactome\\Ensembl2Reactome_All_Levels_2020_07_18.txt");
		final File geneOrderFile = new File("D:\\UMCG\\Genetica\\Projects\\GeneNetwork\\ensgHgncV98.txt");
		final File outputFile = new File("D:\\UMCG\\Genetica\\Projects\\GeneNetwork\\Reactome\\" + pathwayFile.getName() + "_matrix.txt.gz");
		final File outputFile2 = new File("D:\\UMCG\\Genetica\\Projects\\GeneNetwork\\Reactome\\" + pathwayFile.getName() + "_genesInPathways.txt");

		HashMap<String, HashSet<String>> pathwayToGenes = readPathwayFile(pathwayFile);

		ArrayList<String> geneOrder = readGenes(geneOrderFile);

		System.out.println("Total genesets: " + pathwayToGenes.size());

		LinkedHashSet<String> geneOrder2 = new LinkedHashSet<>(geneOrder);

		System.out.println("Genes in order file: " + geneOrder2.size());


		DoubleMatrixDataset<String, String> pathwayMatrix = new DoubleMatrixDataset(geneOrder2, pathwayToGenes.keySet());


		HashSet<String> genesWithPathway = new HashSet<>(10000);
		BufferedWriter geneWriter = new BufferedWriter(new FileWriter(outputFile2));

		for (Map.Entry<String, HashSet<String>> pathwayToGenesEntry : pathwayToGenes.entrySet()) {

			String pathway = pathwayToGenesEntry.getKey();

			for (String gene : pathwayToGenesEntry.getValue()) {

				if (pathwayMatrix.containsRow(gene)) {

					if (genesWithPathway.add(gene)) {
						//add to genes file if not already done
						geneWriter.write(gene);
						geneWriter.write('\n');
					}

					pathwayMatrix.setElement(gene, pathway, 1);
				}

			}

		}

		geneWriter.close();
		pathwayMatrix.save(outputFile);

		System.out.println("Genes in pathway: " + genesWithPathway.size());

	}

	private static HashMap<String, HashSet<String>> readPathwayFile(File pathwayFile) throws Exception {

		final CSVParser parser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();

		CSVReader reader = null;
		if (pathwayFile.getName().endsWith(".gz")) {
			reader = new CSVReaderBuilder(new BufferedReader(
					new InputStreamReader((new GZIPInputStream(new FileInputStream(pathwayFile))))
			)).withSkipLines(0).withCSVParser(parser).build();
		} else {
			reader = new CSVReaderBuilder(new BufferedReader(new FileReader(pathwayFile))).withSkipLines(0).withCSVParser(parser).build();
		}

		HashMap<String, HashSet<String>> pathwayToGenes = new HashMap<>();

		String[] nextLine;
		while ((nextLine = reader.readNext()) != null) {

			if (nextLine[5].equals("Homo sapiens")) {

				String pathway = nextLine[1];
				String ensgId = nextLine[0];

				HashSet<String> pathwayGenes = pathwayToGenes.get(pathway);
				if (pathwayGenes == null) {
					pathwayGenes = new HashSet<>();
					pathwayToGenes.put(pathway, pathwayGenes);
				}

				pathwayGenes.add(ensgId);

			}
		}

		return pathwayToGenes;

	}

	private static ArrayList<String> readGenes(File geneOrderFile) throws IOException {

		final CSVParser parser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
		CSVReader reader = null;
		if (geneOrderFile.getName().endsWith(".gz")) {
			reader = new CSVReaderBuilder(new BufferedReader(new InputStreamReader((new GZIPInputStream(new FileInputStream(geneOrderFile)))))).withSkipLines(0).withCSVParser(parser).build();
		} else {
			reader = new CSVReaderBuilder(new BufferedReader(new FileReader(geneOrderFile))).withSkipLines(0).withCSVParser(parser).build();
		}
		String[] nextLine;
		ArrayList<String> geneOrder = new ArrayList<>();

		while ((nextLine = reader.readNext()) != null) {

			geneOrder.add(nextLine[0]);

		}

		return geneOrder;

	}

}
