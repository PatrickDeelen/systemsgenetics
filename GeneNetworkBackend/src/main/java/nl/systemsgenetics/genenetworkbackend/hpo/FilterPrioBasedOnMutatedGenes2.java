/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.genenetworkbackend.hpo;

import com.opencsv.CSVParser;
import com.opencsv.CSVParserBuilder;
import com.opencsv.CSVReader;
import com.opencsv.CSVReaderBuilder;
import com.opencsv.CSVWriter;

import java.io.*;
import java.util.HashSet;
import java.util.zip.GZIPInputStream;

/**
 * @author patri
 */
public class FilterPrioBasedOnMutatedGenes2 {

	/**
	 * @param args the command line arguments
	 */
	public static void main(String[] args) throws FileNotFoundException, IOException {

//		final File sampleFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\BenchmarkSamples\\Prioritisations\\samplesWithGeno.txt");
//		final File genoFolder = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\BenchmarkSamples\\Prioritisations\\gavinRes\\");
//		final File prioFolder = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\BenchmarkSamples\\Prioritisations");
//		final File resultFolder = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\BenchmarkSamples\\Prioritisations\\rankingCandidateGenes");


//		final File sampleFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\PrioritizeRequests\\Prioritisations\\samples.txt");
//		final File genoFolder = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\PrioritizeRequests\\CandidateGenes\\");
//		final File prioFolder = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\PrioritizeRequests\\Prioritisations");
//		final File resultFolder = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\PrioritizeRequests\\rankingCandidateGenes");

		final File sampleFile = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\BenchmarkSamples\\New5gpm\\hpo5gpm.txt");
		final File genoFolder = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\BenchmarkSamples\\New5gpm\\Genes\\");
		final File prioFolder = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\BenchmarkSamples\\New5gpm\\Prioritisations\\");
		final File resultFolder = new File("C:\\UMCG\\Genetica\\Projects\\GeneNetwork\\BenchmarkSamples\\New5gpm\\RankingCandidateGenes\\");

		final CSVParser parser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
		CSVReader sampleFileReader = null;
		if (sampleFile.getName().endsWith(".gz")) {
			sampleFileReader = new CSVReaderBuilder(new BufferedReader(new InputStreamReader((new GZIPInputStream(new FileInputStream(sampleFile)))))).withSkipLines(0).withCSVParser(parser).build();
		} else {
			sampleFileReader = new CSVReaderBuilder(new BufferedReader(new FileReader(sampleFile))).withSkipLines(0).withCSVParser(parser).build();
		}
		String[] nextLine;
		while ((nextLine = sampleFileReader.readNext()) != null) {

			String sample = nextLine[0];

			String genoSampleName = sample + ".txt";

			File genoFile = new File(genoFolder, genoSampleName);
			File prioFile = new File(prioFolder, sample + ".txt");
			File rankingFile = new File(resultFolder, sample + ".txt");

			System.out.println("------------------------------------------------------------------");
			System.out.println("Sample: " + sample);
			System.out.println("Geno: " + genoFile.getAbsolutePath());
			System.out.println("Prio: " + prioFile.getAbsolutePath());
			System.out.println("Ranking: " + rankingFile.getAbsolutePath());

			HashSet<String> genesWithMutation = getMutatedGenes(genoFile, 0, 0);

			CSVReader prioFileReader = null;
			if (prioFile.getName().endsWith(".gz")) {
				prioFileReader = new CSVReaderBuilder(new BufferedReader(new InputStreamReader((new GZIPInputStream(new FileInputStream(prioFile)))))).withSkipLines(0).withCSVParser(parser).build();
			} else {
				prioFileReader = new CSVReaderBuilder(new BufferedReader(new FileReader(prioFile))).withSkipLines(0).withCSVParser(parser).build();
			}


			CSVWriter writer = new CSVWriter(new FileWriter(rankingFile), '\t', '\0', '\0', "\n");

			String[] outputLine = prioFileReader.readNext();
			writer.writeNext(outputLine);

			while ((outputLine = prioFileReader.readNext()) != null) {

				if (genesWithMutation.contains(outputLine[1])) {
					writer.writeNext(outputLine);
				}

			}

			writer.close();
			prioFileReader.close();

		}

	}

	private static HashSet<String> getMutatedGenes(File genoFile, int colWithGene, int skipHeaderLines) throws IOException {

		final CSVParser parser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
		CSVReader reader = null;
		if (genoFile.getName().endsWith(".gz")) {
			reader = new CSVReaderBuilder(new BufferedReader(new InputStreamReader((new GZIPInputStream(new FileInputStream(genoFile)))))).withSkipLines(skipHeaderLines).withCSVParser(parser).build();
		} else {
			reader = new CSVReaderBuilder(new BufferedReader(new FileReader(genoFile))).withSkipLines(skipHeaderLines).withCSVParser(parser).build();
		}
		HashSet<String> genes = new HashSet<>();

		String[] nextLine;
		while ((nextLine = reader.readNext()) != null) {

			genes.add(nextLine[colWithGene]);

		}

		reader.close();

		return genes;

	}

}
