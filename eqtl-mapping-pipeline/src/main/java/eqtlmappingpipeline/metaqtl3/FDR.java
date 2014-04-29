/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.metaqtl3;

import gnu.trove.map.hash.TDoubleIntHashMap;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.eQTLTextFile;
import umcg.genetica.text.Strings;

/**
 *
 * @author harmjan
 */
public class FDR {

//    public static String permutationDir = null;
//    public static String outputDir = null;
    public enum FDRMethod {

        PROBELEVEL, GENELEVEL, FULL
    };

    public enum FileFormat {

        LARGE, REDUCED
    };

    /**
     * calculate the FalseDiscoveryRate for the discovered eQTLS
     *
     * @param eQTLTextFileLoc the location where the eQTL text files are stored
     * @param nrPermutationsFDR number of permutations performed
     * @param maxNrMostSignificantEQTLs maximum number of eQTLs to output
     * @param fdrcutoff the FDR cutoff
     * @param createQQPlot create a QQ plot after performing FDR calculations
     * @param outputDir set an alternate directory for output
     * @param permutationDir set an alternate directory for permutation files
     * @throws IOException
     */
    public static void calculateFDR(String eQTLTextFileLoc, int nrPermutationsFDR, int maxNrMostSignificantEQTLs, double fdrcutoff, boolean createQQPlot, String outputDir, String permutationDir) throws IOException {

        if (eQTLTextFileLoc == null || eQTLTextFileLoc.length() == 0) {
            throw new IllegalArgumentException("File containing real effects is not specified.");
        }
        if (nrPermutationsFDR < 1) {
            throw new IllegalArgumentException("Need at least one permutation to determine FDR");
        }
        if (maxNrMostSignificantEQTLs < 1) {
            throw new IllegalArgumentException("Need at least a single effect to perform FDR estimation");
        }
        if (fdrcutoff < 0 || fdrcutoff > 1) {
            throw new IllegalArgumentException("FDR threshold should be between 0.0 and 1.0! (Specified: " + fdrcutoff + ")");
        }

        //Load permuted data:
//        // load values for each permutation round:
        if (permutationDir == null) {
            permutationDir = eQTLTextFileLoc;
        }

        if (outputDir == null) {
            outputDir = eQTLTextFileLoc;
        }

        String fileString = permutationDir + "/PermutedEQTLsPermutationRound" + 1 + ".txt.gz";
        TextFile tf = new TextFile(fileString, TextFile.R);
        tf.readLine();
        String[] elems = tf.readLineElems(TextFile.tab);
        int nrColsInPermutedFiles = elems.length;
        tf.close();

        System.out.println(nrColsInPermutedFiles + " columns in permuted QTL file.");
        // new permutationfile format requires different column layout...
        if (nrColsInPermutedFiles > 7) {
            System.out.println("Large permutation files detected.");
            runFDR(eQTLTextFileLoc, nrPermutationsFDR, maxNrMostSignificantEQTLs, fdrcutoff, FileFormat.LARGE, FDRMethod.FULL, outputDir, permutationDir, createQQPlot);
            runFDR(eQTLTextFileLoc, nrPermutationsFDR, maxNrMostSignificantEQTLs, fdrcutoff, FileFormat.LARGE, FDRMethod.PROBELEVEL, outputDir, permutationDir, createQQPlot);
            runFDR(eQTLTextFileLoc, nrPermutationsFDR, maxNrMostSignificantEQTLs, fdrcutoff, FileFormat.LARGE, FDRMethod.GENELEVEL, outputDir, permutationDir, createQQPlot);
        } else {
            runFDR(eQTLTextFileLoc, nrPermutationsFDR, maxNrMostSignificantEQTLs, fdrcutoff, FileFormat.REDUCED, FDRMethod.FULL, outputDir, permutationDir, createQQPlot);
            runFDR(eQTLTextFileLoc, nrPermutationsFDR, maxNrMostSignificantEQTLs, fdrcutoff, FileFormat.REDUCED, FDRMethod.PROBELEVEL, outputDir, permutationDir, createQQPlot);
            if (nrColsInPermutedFiles >= 4) {
                runFDR(eQTLTextFileLoc, nrPermutationsFDR, maxNrMostSignificantEQTLs, fdrcutoff, FileFormat.REDUCED, FDRMethod.GENELEVEL, outputDir, permutationDir, createQQPlot);
            }
        }
    }

    private static void runFDR(String baseDir, int nrPermutationsFDR, int maxNrMostSignificantEQTLs,
            double fdrcutoff, FileFormat f, FDRMethod m, String outputDir, String permutationDir, boolean createQQPlot) throws IOException {
        //Load permuted data:
        // load values for each permutation round:
        System.out.println("");
        if (m == FDRMethod.GENELEVEL) {
            System.out.println("Performing gene level FDR");
        } else if (m == FDRMethod.PROBELEVEL) {
            System.out.println("Performing probe level FDR");
        } else if (m == FDRMethod.FULL) {
            System.out.println("Determining the FDR using all data");
        }


        TDoubleIntHashMap permutedPvalues = new TDoubleIntHashMap(10000, 0.5f);

//        ProgressBar pb = new ProgressBar(nrPermutationsFDR, "Reading permuted data:");
        System.out.println("Reading permuted files");

        for (int permutationRound = 0; permutationRound < nrPermutationsFDR; permutationRound++) {
            String fileString = permutationDir + "/PermutedEQTLsPermutationRound" + (permutationRound + 1) + ".txt.gz";
            System.out.println(fileString);
            // read the permuted eqtl output
            TextFile gz = new TextFile(fileString, TextFile.R);


            String[] header = gz.readLineElems(TextFile.tab);
            int snpcol = -1;
            int pvalcol = -1;
            int probecol = -1;
            int genecol = -1;

            if (f == FileFormat.REDUCED) {



                //PValue  SNP     Probe   Gene 
                for (int col = 0; col < header.length; col++) {
                    if (header[col].equals("PValue")) {
                        pvalcol = col;
                    }
                    if (header[col].equals("SNP")) {
                        snpcol = col;
                    }
                    if (header[col].equals("Probe")) {
                        probecol = col;
                    }
                    if (header[col].equals("Gene")) {
                        genecol = col;
                    }
                }

                //PValue  SNP     Probe   Gene
                if (snpcol == -1 || pvalcol == -1 || probecol == -1 && genecol == -1) {
                    System.out.println("Column not found in permutation file: " + fileString);
                    System.out.println("PValue: " + pvalcol);
                    System.out.println("SNP: " + snpcol);
                    System.out.println("Probe: " + probecol);
                    System.out.println("Gene: " + genecol);
                }
            }
            String[] data = gz.readLineElemsReturnReference(TextFile.tab);
            int itr = 0;

            HashSet<String> visitedEffects = new HashSet<String>();
            while (data != null) {

                if (data.length != 0) {
                    if (itr > maxNrMostSignificantEQTLs - 1) {
                        System.out.println("Breaking because: " + itr);
                        break;
                    } else {
                        int filteronColumn;
                        String fdrId = null;
                        if (f == FileFormat.REDUCED) {
                            if (m == FDRMethod.FULL) {
                                //fdrId = data[snpcol] + "-" + data[probecol];
                                filteronColumn = probecol;
                            } else if (m == FDRMethod.GENELEVEL && data.length > 3) {
                                fdrId = data[genecol];
                                filteronColumn = genecol;
                            } else {
                                fdrId = data[probecol];
                                filteronColumn = probecol;
                            }

                        } else {
                            if (m == FDRMethod.GENELEVEL) {
                                fdrId = data[eQTLTextFile.HUGO];
                                filteronColumn = eQTLTextFile.HUGO;
                            } else if (m == FDRMethod.PROBELEVEL) {
                                fdrId = data[4];
                                filteronColumn = 4;
                            }
                        }

                        // take top effect per gene / probe
                        if (m == FDRMethod.FULL || (!fdrId.equals("-") && !visitedEffects.contains(fdrId))) {

                            if (m != FDRMethod.FULL) {
                                visitedEffects.add(fdrId);
                            }

                            double permutedP = Double.parseDouble(data[0]);
                            if (permutedPvalues.containsKey(permutedP)) {
                                permutedPvalues.increment(permutedP);
                            } else {
                                permutedPvalues.put(permutedP, 1);
                            }

                            itr++;
                        }

                        data = gz.readLineElemsReturnReference(TextFile.tab);
                    }
                }
            }
            gz.close();


        }




        double[] uniquePermutedPvalues = permutedPvalues.keys();
        Arrays.sort(uniquePermutedPvalues);

        double[] uniquePermutedPvaluesCounts = new double[uniquePermutedPvalues.length];

        long cummulativeCount = 0;
        double nrPermutationsFDRd = (double) nrPermutationsFDR;
        for (int i = 0; i < uniquePermutedPvalues.length; ++i) {

            cummulativeCount += permutedPvalues.get(uniquePermutedPvalues[i]);
            uniquePermutedPvaluesCounts[i] = cummulativeCount / nrPermutationsFDRd;

        }
        permutedPvalues = null;
        System.out.println("Number of unique permutation p-values: " + uniquePermutedPvalues.length);


        String outFileName;
        String outFileNameSnps;
        String outFileNameProbes;
        String outFileNameAll;

        if (outputDir == null) {
            outputDir = baseDir;
        }

        if (m == FDRMethod.GENELEVEL) {
            outFileName = outputDir + "/eQTLsFDR" + fdrcutoff + "-GeneLevel.txt";
            outFileNameSnps = outputDir + "/eQTLSNPsFDR" + fdrcutoff + "-GeneLevel.txt";
            outFileNameProbes = outputDir + "/eQTLProbesFDR" + fdrcutoff + "-GeneLevel.txt";
            outFileNameAll = outputDir + "/eQTLsFDR-GeneLevel.txt.gz";
        } else if (m == FDRMethod.PROBELEVEL) {
            outFileName = outputDir + "/eQTLsFDR" + fdrcutoff + "-ProbeLevel.txt";
            outFileNameSnps = outputDir + "/eQTLSNPsFDR" + fdrcutoff + "-ProbeLevel.txt";
            outFileNameProbes = outputDir + "/eQTLProbesFDR" + fdrcutoff + "-ProbeLevel.txt";
            outFileNameAll = outputDir + "/eQTLsFDR-ProbeLevel.txt.gz";
        } else {
            outFileName = outputDir + "/eQTLsFDR" + fdrcutoff + ".txt";
            outFileNameSnps = outputDir + "/eQTLSNPsFDR" + fdrcutoff + ".txt";
            outFileNameProbes = outputDir + "/eQTLProbesFDR" + fdrcutoff + ".txt";
            outFileNameAll = outputDir + "/eQTLsFDR.txt.gz";
        }

        BufferedWriter outputWriterSignificant = new BufferedWriter(new FileWriter(outFileName));
        BufferedWriter outputWriterESNPs = new BufferedWriter(new FileWriter(outFileNameSnps));
        BufferedWriter outputWriterEProbes = new BufferedWriter(new FileWriter(outFileNameProbes));
        BufferedWriter outputWriterAll = new BufferedWriter(new FileWriter(outFileNameAll));

        String fileString = baseDir + "/eQTLs.txt.gz";
        if (!Gpio.exists(fileString)) {
            System.out.println("Could not find file: " + fileString + " trying un-GZipped file....");
            fileString = baseDir + "/eQTLs.txt";
        }
        if (!Gpio.exists(fileString)) {
            System.out.println("Could not find file: " + fileString);
            System.exit(0);
        }


        TextFile realEQTLs = new TextFile(fileString, TextFile.R);

        String header = realEQTLs.readLine();

        outputWriterAll.append(header);
        outputWriterAll.append("\tFDR\n");

        outputWriterEProbes.append(header);
        outputWriterEProbes.append("\tFDR\n");
        
        outputWriterESNPs.append(header);
        outputWriterESNPs.append("\tFDR\n");
        
        outputWriterSignificant.append(header);
        outputWriterSignificant.append("\tFDR\n");

        String str = realEQTLs.readLine();

// REAL DATA PROCESSING
        int itr = 0;
        HashSet<String> visitedEffects = new HashSet<String>();
        HashSet<String> visitedSnps = new HashSet<String>();
        HashSet<String> visitedProbes = new HashSet<String>();
        double lastEqtlPvalue = 0;

        double currentPvalue = 0;
        ArrayList<String> currentPvalueEqtls = new ArrayList<String>();
        ArrayList<String> currentPvalueEqtlSnps = new ArrayList<String>();
        ArrayList<String> currentPvalueEqtlProbes = new ArrayList<String>();


        int lastUsedPermutedPvalueIndex = 0;


        int nrSignificantEQTLs = 0;

        while (str != null) {
            if (itr > maxNrMostSignificantEQTLs - 1) {
                break;
            } else {

                String fdrId = null;
                String[] data = Strings.tab.split(str);

                if (m == FDRMethod.GENELEVEL) {
                    fdrId = data[eQTLTextFile.HUGO];
                } else if (m == FDRMethod.PROBELEVEL) {
                    fdrId = data[4];
                }

                double eQtlPvalue = Double.parseDouble(data[0]);

                if (itr > 0 && lastEqtlPvalue > eQtlPvalue) {
                    System.err.println("Sorted P-Value list is not perfectly sorted!!!!");
                    System.exit(-1);
                }

                if (eQtlPvalue > currentPvalue) {
                    //Process old results for current pvalue

                    double fdr = 0;
                    if (currentPvalue >= uniquePermutedPvalues[0]) {

                        while (uniquePermutedPvalues[lastUsedPermutedPvalueIndex + 1] <= currentPvalue && lastUsedPermutedPvalueIndex < uniquePermutedPvalues.length - 2) {
                            ++lastUsedPermutedPvalueIndex;
                        }
                        fdr = uniquePermutedPvaluesCounts[lastUsedPermutedPvalueIndex] / itr;

                    }

                    for (int i = 0; i < currentPvalueEqtls.size(); ++i) {
                        String cachedEqtls = currentPvalueEqtls.get(i);
                        String cachedEqtlsProbe = currentPvalueEqtlProbes.get(i);
                        String cachedEqtlsSnps = currentPvalueEqtlSnps.get(i);

                        StringBuilder currentString = new StringBuilder();
                        currentString.append(cachedEqtls).append('\t').append(String.valueOf(fdr)).append('\n');

                        outputWriterAll.append(currentString.toString());

                        if (fdr <= fdrcutoff) {
                            if (!visitedProbes.contains(cachedEqtlsProbe)) {
                                outputWriterEProbes.append(currentString.toString());
                                visitedProbes.add(cachedEqtlsProbe);
                            }
                            if (!visitedSnps.contains(cachedEqtlsSnps)) {
                                outputWriterESNPs.append(currentString.toString());
                                visitedSnps.add(cachedEqtlsSnps);

                            }

                            outputWriterSignificant.append(currentString.toString());
                            ++nrSignificantEQTLs;
                        }

                    }


                    //Create new temp list for this pvalue
                    currentPvalue = eQtlPvalue;
                    currentPvalueEqtls.clear();
                    currentPvalueEqtlProbes.clear();
                    currentPvalueEqtlSnps.clear();
                    currentPvalueEqtls.add(str);
                    currentPvalueEqtlProbes.add(data[eQTLTextFile.PROBE]);
                    currentPvalueEqtlSnps.add(data[eQTLTextFile.SNP]);

                } else {
                    //add to current pvalue list
                    currentPvalueEqtls.add(str);
                    currentPvalueEqtlProbes.add(data[eQTLTextFile.PROBE]);
                    currentPvalueEqtlSnps.add(data[eQTLTextFile.SNP]);
                }

                lastEqtlPvalue = eQtlPvalue;
                
                if (m == FDRMethod.FULL || (!fdrId.equals("-") && !visitedEffects.contains(fdrId))) {
                    itr++;
                    visitedEffects.add(fdrId);
                }

                str = realEQTLs.readLine();
            }

        }

        //Write buffer to files
        double fdr = 0;
        if (currentPvalue >= uniquePermutedPvalues[0]) {

            while (uniquePermutedPvalues[lastUsedPermutedPvalueIndex + 1] <= currentPvalue && lastUsedPermutedPvalueIndex < uniquePermutedPvalues.length - 2) {
                ++lastUsedPermutedPvalueIndex;
            }
            fdr = uniquePermutedPvaluesCounts[lastUsedPermutedPvalueIndex] / itr;

        }

        for (int i = 0; i < currentPvalueEqtls.size(); ++i) {
            String cachedEqtls = currentPvalueEqtls.get(i);
            String cachedEqtlsProbe = currentPvalueEqtlProbes.get(i);
            String cachedEqtlsSnps = currentPvalueEqtlSnps.get(i);

            StringBuilder currentString = new StringBuilder();
            currentString.append(cachedEqtls).append('\t').append(String.valueOf(fdr)).append('\n');

            outputWriterAll.append(currentString.toString());

            if (fdr <= fdrcutoff) {
                if (!visitedProbes.contains(cachedEqtlsProbe)) {
                    outputWriterEProbes.append(currentString.toString());
                    visitedSnps.add(cachedEqtlsProbe);
                }
                if (!visitedSnps.contains(cachedEqtlsSnps)) {
                    outputWriterESNPs.append(currentString.toString());
                    visitedSnps.add(cachedEqtlsSnps);

                }

                outputWriterSignificant.append(currentString.toString());
                ++nrSignificantEQTLs;
            }

        }


        realEQTLs.close();
        outputWriterAll.close();
        outputWriterEProbes.close();
        outputWriterESNPs.close();
        outputWriterSignificant.close();

        //System.out.println("");
        System.out.println("Number of significant eQTLs:\t" + nrSignificantEQTLs);
        System.out.println(" - Number of unique SNPs, constituting an eQTL:\t" + visitedSnps.size());
        System.out.println(" - Number of unique probes, constituting an eQTL:\t" + visitedProbes.size());

        if (createQQPlot) {

            System.err.println("Sorry, QQ plot function is temporarily (or for a very long time) unavailable.");

            //System.out.println("Creating QQ plot. This might take a while...");
            //QQPlot qq = new QQPlot();
            //String fileName = baseDir + "/eQTLsFDR" + fdrcutoff + fileSuffix + "-QQPlot.pdf";
            //qq.draw(fileName, fdrcutoff, nrPermutationsFDR,
            //		maxNrMostSignificantEQTLs, permutedPValues.toArray(), nrRealDataEQTLs, pValues,
            //		pValueSignificant, nrSignificantEQTLs);
        }

    }
}
