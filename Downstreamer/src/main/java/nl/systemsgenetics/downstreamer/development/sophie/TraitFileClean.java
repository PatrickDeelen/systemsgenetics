/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.downstreamer.development.sophie;

/**
 *
 * @author Sophie Mulc
 *
 */
import com.opencsv.CSVParser;
import com.opencsv.CSVParserBuilder;
import com.opencsv.CSVReader;
import com.opencsv.CSVReaderBuilder;
import com.opencsv.CSVWriter;

import java.io.*;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.zip.GZIPInputStream;

@Deprecated
public class TraitFileClean {

    /**
     * @param args the command line arguments
     * @throws java.io.FileNotFoundException
     */
//        MAKE TRAIT FILE
    static void trait(List<String> iids, File traitFile) throws FileNotFoundException, IOException {

        CSVWriter writer = new CSVWriter(new FileWriter(traitFile), '\t', '\0', '\0', "\n");

        String[] outLine = new String[iids.size() + 1];
        int c = 0;
        outLine[c++] = "PHENO";
        for (String iid : iids) {
            outLine[c++] = iid;
        }
        writer.writeNext(outLine);

        Random randomno1 = new Random();

        for (int j = 1; j < 1001; ++j) {
            c = 0;

            outLine[c++] = "PH" + j;
            for (int i = 0; i < iids.size(); ++i) {
                outLine[c++] = String.valueOf(randomno1.nextGaussian());
            }
            writer.writeNext(outLine);
        }
        writer.close();
    }

//  MAKE PROBE ANNOTATION FILE   
    static void probeAnnotation(File probeAnnotationFile) throws IOException {

        CSVWriter writer1 = new CSVWriter(new FileWriter(probeAnnotationFile), '\t', '\0', '\0', "\n");

        String[] output1Line = new String[6];
        int c = 0;
        output1Line[c++] = "Platform";
        output1Line[c++] = "HT12v4-ArrayAddress";
        output1Line[c++] = "Symbol";
        output1Line[c++] = "Chr";
        output1Line[c++] = "ChrStart";
        output1Line[c++] = "ChrEnd";
        writer1.writeNext(output1Line);

        for (int j = 1; j < 1001; ++j) {
            c = 0;

            output1Line[c++] = "Plat";
            output1Line[c++] = "PH" + j;
            output1Line[c++] = "PH" + j;
            output1Line[c++] = "25";
            output1Line[c++] = "100";
            output1Line[c++] = "154";
            writer1.writeNext(output1Line);
        }
        writer1.close();
    }

    //  MAKE COUPLING FILE    
    static void coupling(List<String> iids, File couplingFile) throws IOException {

        CSVWriter writer2 = new CSVWriter(new FileWriter(couplingFile), '\t', '\0', '\0', "\n");

        String[] output2Line = new String[2];
        for (String iid : iids) {
            int c = 0;
            output2Line[c++] = iid;
            output2Line[c++] = iid;
            writer2.writeNext(output2Line);
        }
        writer2.close();
    }

    public static void main(String[] args) throws FileNotFoundException, IOException {
        // TODO code application logic here 
        File phase3File = new File("C:\\Users\\Sophie Mulc\\Documents\\DEPICT2\\phase3_corrected.psam");
        File traitFile = new File("C:\\Users\\Sophie Mulc\\Documents\\DEPICT2\\TraitFile.txt");
        File probeAnnotationFile = new File("C:\\Users\\Sophie Mulc\\Documents\\DEPICT2\\ProbeAnnotationFile.txt");
        File couplingFile = new File("C:\\Users\\Sophie Mulc\\Documents\\DEPICT2\\CouplingFile.txt");
        //FileReader(String phase3_corrected)
        final CSVParser gmtParser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
        CSVReader gmtReader = null;
        if (phase3File.getName().endsWith(".gz")) {
            gmtReader = new CSVReaderBuilder((new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(phase3File)))))).withSkipLines(1).withCSVParser(gmtParser).build();
        } else {
            gmtReader = new CSVReaderBuilder(new BufferedReader(new FileReader(phase3File))).withSkipLines(1).withCSVParser(gmtParser).build();
        }

        List<String> iids = new ArrayList<>();

        String[] inputLine;
        while ((inputLine = gmtReader.readNext()) != null) {

            String iid = inputLine[0];

            iids.add(iid);
        }

        trait(iids, traitFile);
        probeAnnotation(probeAnnotationFile);
        coupling(iids, couplingFile);
    }

}
