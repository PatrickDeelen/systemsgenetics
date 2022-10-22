package nl.systemsgenetics.genenetworkbackend;

import com.opencsv.CSVParser;
import com.opencsv.CSVParserBuilder;
import com.opencsv.CSVReader;
import com.opencsv.CSVReaderBuilder;

import java.io.*;
import java.util.ArrayList;
import java.util.EnumMap;
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
public class ConvertGoToMatrix {

	private static enum GoType {
		F, P, C
	}

	;

	/**
	 * @param args the command line arguments
	 * @throws java.io.IOException
	 * @throws java.lang.Exception
	 */
	public static void main(String[] args) throws IOException, Exception {

		final File pathwayFile = new File("D:\\UMCG\\Genetica\\Projects\\GeneNetwork\\GO\\goa_human_2020_06_01.gaf");
		final File geneOrderFile = new File("D:\\UMCG\\Genetica\\Projects\\GeneNetwork\\ensgHgncV98.txt");
		final File uniPortToEnsgMapFile = new File("D:\\UMCG\\Genetica\\Projects\\GeneNetwork\\ensgUniProtIdV98.txt");
		final File hgncToEnsgMapFile = new File("D:\\UMCG\\Genetica\\Projects\\GeneNetwork\\ensgHgncV98.txt");
		final File outputFolder = new File("D:\\UMCG\\Genetica\\Projects\\GeneNetwork\\GO");

		HashMap<String, ArrayList<String>> uniProtToEnsgMap = loadUniProtToEnsgMap(uniPortToEnsgMapFile);
		HashMap<String, ArrayList<String>> hgncToEnsgMap = loadHgncToEnsgMap(hgncToEnsgMapFile);

		EnumMap<GoType, HashMap<String, HashSet<String>>> goTypePathwayToGenes = readPathwayFile(pathwayFile, uniProtToEnsgMap, hgncToEnsgMap);

		ArrayList<String> geneOrder = readGenes(geneOrderFile);

		LinkedHashSet<String> geneOrder2 = new LinkedHashSet<>(geneOrder);

		System.out.println("Genes in order file: " + geneOrder2.size());

		for (GoType goType : GoType.values()) {

			final File outputFile = new File(outputFolder, pathwayFile.getName() + "_" + goType.toString() + "_2020_06_01_matrix.txt.gz");
			final File outputFile2 = new File(outputFolder, pathwayFile.getName() + "_" + goType.toString() + "_2020_06_01_genesInPathways.txt");

			HashMap<String, HashSet<String>> pathwayToGenes = goTypePathwayToGenes.get(goType);

			System.out.println("Total genesets of " + goType.toString() + ": " + pathwayToGenes.size());

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

			pathwayMatrix.save(outputFile);
			geneWriter.close();

			System.out.println("Genes in " + goType.toString() + " pathway: " + genesWithPathway.size());

		}

	}

	private static EnumMap<GoType, HashMap<String, HashSet<String>>> readPathwayFile(File pathwayFile, HashMap<String, ArrayList<String>> uniProtToEnsgMap, HashMap<String, ArrayList<String>> hgncToEnsgMap) throws Exception {

		final CSVParser parser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
		CSVReader reader = null;
		if (pathwayFile.getName().endsWith(".gz")) {
			reader = new CSVReaderBuilder(new BufferedReader(
					new InputStreamReader((new GZIPInputStream(new FileInputStream(pathwayFile))))
			)).withSkipLines(0).withCSVParser(parser).build();
		} else {
			reader = new CSVReaderBuilder(new BufferedReader(new FileReader(pathwayFile))).withSkipLines(0).withCSVParser(parser).build();
		}

		EnumMap<GoType, HashMap<String, HashSet<String>>> goTypePathwayToGenes = new EnumMap<>(GoType.class);

		for (GoType goType : GoType.values()) {
			goTypePathwayToGenes.put(goType, new HashMap<>());
		}

		String[] nextLine;
		while ((nextLine = reader.readNext()) != null) {

			if (nextLine[0].charAt(0) == '!') {
				continue;
			}

			//Exclude special tuples like empty
			if (nextLine[3].isEmpty()) {

				String uniProtId = nextLine[1];
				String hgcnId = nextLine[2];
				String pathway = nextLine[4];

				ArrayList<String> ensgIds = uniProtToEnsgMap.get(uniProtId);

				if (ensgIds == null) {
					ensgIds = hgncToEnsgMap.get(hgcnId);
				}

				if (ensgIds == null) {
					System.err.println("Missing mapping for gene: " + uniProtId + " " + nextLine[2]);
				} else {

					HashMap<String, HashSet<String>> pathwayToGenes = goTypePathwayToGenes.get(GoType.valueOf(nextLine[8]));

					HashSet<String> pathwayGenes = pathwayToGenes.get(pathway);
					if (pathwayGenes == null) {
						pathwayGenes = new HashSet<>();
						pathwayToGenes.put(pathway, pathwayGenes);
					}

					for (String ensgId : ensgIds) {
						pathwayGenes.add(ensgId);
					}

				}

			}

		}

		return goTypePathwayToGenes;

	}

	public static ArrayList<String> readGenes(File geneOrderFile) throws IOException {

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

	private static HashMap<String, ArrayList<String>> loadUniProtToEnsgMap(File map) throws FileNotFoundException, IOException, Exception {

		final CSVParser parser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
		CSVReader reader = null;
		if (map.getName().endsWith(".gz")) {
			reader = new CSVReaderBuilder(new BufferedReader(
					new InputStreamReader((new GZIPInputStream(new FileInputStream(map))))
			)).withSkipLines(0).withCSVParser(parser).build();
		} else {
			reader = new CSVReaderBuilder(new BufferedReader(new FileReader(map))).withSkipLines(0).withCSVParser(parser).build();
		}
		String[] nextLine = reader.readNext();

//		if (!nextLine[0].equals("Gene stable ID") || !nextLine[1].equals("Transcript stable ID") || !nextLine[2].equals("UniProtKB Gene Name ID") || !nextLine[3].equals("UniProtKB/Swiss-Prot ID") || !nextLine[4].equals("UniProtKB/TrEMBL ID")) {
//			throw new Exception("Header of uniprot to ensg map should be: \"Gene stable ID[tab]Transcript stable ID[tab]UniProtKB Gene Name ID[tab]UniProtKB/Swiss-Prot ID[tab]UniProtKB/TrEMBL ID\"");
//		}
		if (!nextLine[0].equals("Gene stable ID") || !nextLine[1].equals("Transcript stable ID") || !nextLine[2].equals("UniProtKB/Swiss-Prot ID")) {
			throw new Exception("Header of uniprot to ensg map should be: \"Gene stable ID[tab]Transcript stable ID[tab]UniProtKB/Swiss-Prot ID\"");
		}

		HashMap<String, ArrayList<String>> uniProtToEnsgMap = new HashMap<>(70000);

		while ((nextLine = reader.readNext()) != null) {

			String uniProtId = nextLine[2];

			ArrayList<String> uniProtEnsgIds = uniProtToEnsgMap.get(uniProtId);
			if (uniProtEnsgIds == null) {
				uniProtEnsgIds = new ArrayList<>();
				uniProtToEnsgMap.put(uniProtId, uniProtEnsgIds);
			}

			uniProtEnsgIds.add(nextLine[0]);

		}

		return uniProtToEnsgMap;

	}

	private static HashMap<String, ArrayList<String>> loadHgncToEnsgMap(File map) throws FileNotFoundException, IOException, Exception {

		final CSVParser parser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
		CSVReader reader = null;
		if (map.getName().endsWith("gz")) {
			reader = new CSVReaderBuilder(new BufferedReader(new InputStreamReader((new GZIPInputStream(new FileInputStream(map)))))).withSkipLines(0).withCSVParser(parser).build();

		} else {
			reader = new CSVReaderBuilder(new BufferedReader(new FileReader(map))).withSkipLines(0).withCSVParser(parser).build();
		}


		String[] nextLine = reader.readNext();

//		if (!nextLine[0].equals("Gene stable ID") || !nextLine[1].equals("Transcript stable ID") || !nextLine[2].equals("UniProtKB Gene Name ID") || !nextLine[3].equals("UniProtKB/Swiss-Prot ID") || !nextLine[4].equals("UniProtKB/TrEMBL ID")) {
//			throw new Exception("Header of uniprot to ensg map should be: \"Gene stable ID[tab]Transcript stable ID[tab]UniProtKB Gene Name ID[tab]UniProtKB/Swiss-Prot ID[tab]UniProtKB/TrEMBL ID\"");
//		}
		if (!nextLine[0].equals("Gene stable ID") || !nextLine[1].equals("HGNC symbol")) {
			throw new Exception("Header of hgnc to ensg map should be: \"Gene stable ID[tab]HGNC symbol\"");
		}

		HashMap<String, ArrayList<String>> hgncToEnsgMap = new HashMap<>(70000);

		while ((nextLine = reader.readNext()) != null) {

			String hgncId = nextLine[1];

			ArrayList<String> hgncEnsgIds = hgncToEnsgMap.get(hgncId);
			if (hgncEnsgIds == null) {
				hgncEnsgIds = new ArrayList<>();
				hgncToEnsgMap.put(hgncId, hgncEnsgIds);
			}

			hgncEnsgIds.add(nextLine[0]);

		}

		return hgncToEnsgMap;

	}

}
