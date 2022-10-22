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
public class ConvertGmtToMatrix {

	/**
	 * @param args the command line arguments
	 * @throws java.io.IOException
	 * @throws java.lang.Exception
	 */
	public static void main(String[] args) throws IOException, Exception {

		final File gmtFile = new File("D:\\UMCG\\Genetica\\Projects\\GeneNetwork\\KEGG\\c2.cp.kegg.v7.1.entrez.gmt");
		final File ncbiToEnsgMapFile = new File("D:\\UMCG\\Genetica\\Projects\\GeneNetwork\\ensgNcbiIdV98.txt");
		final File geneOrderFile = new File("D:\\UMCG\\Genetica\\Projects\\GeneNetwork\\ensgNcbiIdV98.txt");
		final File outputFile = new File("D:\\UMCG\\Genetica\\Projects\\GeneNetwork\\KEGG\\" + gmtFile.getName() + "_matrix.txt.gz");
		final File outputFile2 = new File("D:\\UMCG\\Genetica\\Projects\\GeneNetwork\\KEGG\\" + gmtFile.getName() + "_genesInPathways.txt");

		HashMap<String, String> ncbiToEnsgMap = loadNcbiToEnsgMap(ncbiToEnsgMapFile);

		HashMap<String, HashSet<String>> gmtPathwayToGenes = readGmtFile(gmtFile, ncbiToEnsgMap);

		ArrayList<String> geneOrder = readGenes(geneOrderFile);
		LinkedHashSet<String> geneOrder2 = new LinkedHashSet<>(geneOrder);

		System.out.println("Total genesets: " + gmtPathwayToGenes.size());
		System.out.println("Genes in order file: " + geneOrder2.size());

		DoubleMatrixDataset<String, String> gmtMatrix = new DoubleMatrixDataset(geneOrder2, gmtPathwayToGenes.keySet());

		HashSet<String> genesWithPathway = new HashSet<>(10000);
		BufferedWriter geneWriter = new BufferedWriter(new FileWriter(outputFile2));

		for (Map.Entry<String, HashSet<String>> gmtPathwayToGenesEntry : gmtPathwayToGenes.entrySet()) {

			String gmtPathway = gmtPathwayToGenesEntry.getKey();

			for (String gene : gmtPathwayToGenesEntry.getValue()) {

				if (gmtMatrix.containsRow(gene)) {

					if (genesWithPathway.add(gene)) {
						//add to genes file if not already done
						geneWriter.write(gene);
						geneWriter.write('\n');
					}

					if (!gmtMatrix.containsCol(gmtPathway)) {
						System.out.println("missing pathway: " + gmtPathway);
					}

					if (!gmtMatrix.containsRow(gene)) {
						System.out.println("Missing gene: " + gene);
					}

					gmtMatrix.setElement(gene, gmtPathway, 1);

				}
			}

		}

		gmtMatrix.save(outputFile);
		geneWriter.close();

		System.out.println("Genes in pathway: " + genesWithPathway.size());

	}

	private static HashMap<String, String> loadNcbiToEnsgMap(File ncbiToEnsgMapFile) throws FileNotFoundException, IOException, Exception {

		final CSVParser parser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
		CSVReader reader = null;
		if (ncbiToEnsgMapFile.getName().endsWith(".gz")) {
			reader = new CSVReaderBuilder(new BufferedReader(new FileReader(ncbiToEnsgMapFile))).withSkipLines(0).withCSVParser(parser).build();
		} else {
			reader = new CSVReaderBuilder(new BufferedReader(new FileReader(ncbiToEnsgMapFile))).withSkipLines(0).withCSVParser(parser).build();
		}
		String[] nextLine = reader.readNext();

		if (!nextLine[0].equals("Gene stable ID") || !nextLine[1].equals("NCBI gene ID")) {
			throw new Exception("Header of ncbi to ensg map should be: \"Gene stable ID	NCBI gene ID\"");
		}

		HashMap<String, String> ncbiToEnsgMap = new HashMap<>(70000);

		while ((nextLine = reader.readNext()) != null) {
			ncbiToEnsgMap.put(nextLine[1], nextLine[0]);
		}

		return ncbiToEnsgMap;

	}

	private static HashMap<String, HashSet<String>> readGmtFile(File hpoFile, HashMap<String, String> ncbiToEnsgMap) throws Exception {

		final CSVParser gmtParser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
		CSVReader gmtReader = null;
		if (hpoFile.getName().endsWith(".gz")) {
			gmtReader = new CSVReaderBuilder(new BufferedReader(
					new InputStreamReader((new GZIPInputStream(new FileInputStream(hpoFile))))
			)).withSkipLines(0).withCSVParser(gmtParser).build();
		} else {
			gmtReader = new CSVReaderBuilder(new BufferedReader(new FileReader(hpoFile))).withSkipLines(0).withCSVParser(gmtParser).build();
		}


		HashMap<String, HashSet<String>> gmtPathwayToGenes = new HashMap<>();

		String[] nextLine;
		while ((nextLine = gmtReader.readNext()) != null) {

			String gmtPathway = nextLine[0];

			if (gmtPathwayToGenes.containsKey(gmtPathway)) {
				throw new Exception("Found pathway twice in GMT file: " + gmtPathway);
			}

			HashSet<String> gmtGenes = new HashSet<>();
			gmtPathwayToGenes.put(gmtPathway, gmtGenes);

			for (int i = 2; i < nextLine.length; ++i) {

				String ncbiId = nextLine[i];
				String ensgId = ncbiToEnsgMap.get(ncbiId);
				if (ensgId == null) {
					System.err.println("Missing mapping for gene: " + ncbiId + " " + nextLine[1]);
				} else {

					gmtGenes.add(ensgId);

				}

			}

		}

		return gmtPathwayToGenes;

	}

	private static ArrayList<String> readGenes(File geneOrderFile) throws IOException {

		final CSVParser parser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
		CSVReader reader = null;
		if (geneOrderFile.getName().endsWith(".gz")) {
			reader = new CSVReaderBuilder(new BufferedReader(
					new InputStreamReader((new GZIPInputStream(new FileInputStream(geneOrderFile))))
			)).withSkipLines(0).withCSVParser(parser).build();
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
