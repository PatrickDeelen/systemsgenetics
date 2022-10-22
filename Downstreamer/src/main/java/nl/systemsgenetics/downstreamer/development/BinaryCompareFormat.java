/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.downstreamer.development;

import com.opencsv.CSVParser;
import com.opencsv.CSVParserBuilder;
import com.opencsv.CSVReader;
import com.opencsv.CSVReaderBuilder;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;

import java.io.*;
import java.util.zip.GZIPInputStream;

/**
 * @author Sophie Mulc
 */
public class BinaryCompareFormat {

	/**
	 * @param args the command line arguments
	 * @throws java.io.IOException
	 */
	public static void main(String[] args) throws IOException {
		// TODO code application logic here

		//To compare the format of the files we created in the eQTL pipeline for running GWAS to the format needed to run in Lude's code.
		final File predictionMatrixFile = new File("/groups/umcg-wijmenga/scr02/umcg-smulcahy/eQTLResults_texttrue2/eQTL.binary");
		System.out.println(predictionMatrixFile);
		DoubleMatrixDataset<String, String> predictionMatrixFull = DoubleMatrixDataset.loadTransEqtlExpressionMatrix(predictionMatrixFile.getAbsolutePath());

		System.out.println(predictionMatrixFull.getElement("rs351365", "PH443"));
		File eQTLfile = new File("/groups/umcg-wijmenga/scr02/umcg-smulcahy/eQTLResults_texttrue2/eQTLs.txt.gz");
		final CSVParser eQTLparser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
		CSVReader eQTLreader = null;
		if (eQTLfile.getName().endsWith(".gz")) {
			eQTLreader = new CSVReaderBuilder((new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(eQTLfile)))))).withSkipLines(1).withCSVParser(eQTLparser).build();
		} else {
			eQTLreader = new CSVReaderBuilder(new BufferedReader(new FileReader(eQTLfile))).withCSVParser(eQTLparser).withSkipLines(1).build();
		}


		int c = 0;
		String[] inputLine;
		while ((inputLine = eQTLreader.readNext()) != null) {

			String snp = inputLine[1];
			String pheno = inputLine[4];
			String zscore = inputLine[10];

			// Test with one site only need this line and to read in file: predictionMatrixFull.getElement("rs351365", "PH443"); i.e (rowName(snp), columnName(pheno))


			double zscore2 = (predictionMatrixFull.getElement(snp, pheno));

			//convert string to double, then look up how to compare two doubles - this is to compare the zscores
			double d_zscore = Double.parseDouble(zscore);

			double compare_zscores = d_zscore - zscore2;

			//count occurances of above 0 comparisons
			if (Math.abs(compare_zscores) > 0.00001) {
				c++;
			}
		}
		System.out.println("Number of occurrances where z-scores differ: " + c);
		//save new file
		predictionMatrixFull.saveBinary("/groups/umcg-wijmenga/scr02/umcg-smulcahy/eQTLResults_texttrue2/eQTL2");

		System.out.println("Done saving");

	}

}
