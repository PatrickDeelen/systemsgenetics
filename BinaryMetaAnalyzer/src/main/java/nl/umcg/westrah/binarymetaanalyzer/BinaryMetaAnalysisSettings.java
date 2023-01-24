/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.westrah.binarymetaanalyzer;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.NoSuchElementException;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.apache.commons.configuration.ConfigurationException;
import org.apache.commons.configuration.XMLConfiguration;
import org.apache.commons.lang3.RandomStringUtils;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;

/**
 * @author harm-jan
 */
public class BinaryMetaAnalysisSettings {


    @Override
    public String toString() {
        return "BinaryMetaAnalysisSettings{" +
                "nrPermutations=" + nrPermutations +
                "\n startPermutations=" + startPermutations +
                "\n useAbsoluteZscore=" + useAbsoluteZscore +
                "\n finalEQTLBufferMaxLength=" + finalEQTLBufferMaxLength +
                "\n nrOfBins=" + nrOfBins +
                "\n includeSNPsWithoutProperMapping=" + includeSNPsWithoutProperMapping +
                "\n includeProbesWithoutProperMapping=" + includeProbesWithoutProperMapping +
                "\n analysisType=" + analysisType +
                "\n cisdistance=" + cisdistance +
                "\n transdistance=" + transdistance +
                "\n makezscoreplot=" + makezscoreplot +
                "\n probetranslationfile='" + probetranslationfile + '\'' +
                "\n rescalingOfSampleSize=" + rescalingOfSampleSize +
                "\n featureOccuranceScaleMaps=" + featureOccuranceScaleMaps +
                "\n datasetnames=" + datasetnames +
                "\n datasetPrefix=" + datasetPrefix +
                "\n datasetlocations=" + datasetlocations +
                "\n datasetannotations=" + datasetannotations +
                "\n selectedProbes=" + selectedProbes +
                "\n output='" + output + '\'' +
                "\n makezscoretable=" + makezscoretable +
                "\n confineSNPs=" + confineSNPs +
                "\n probeDatasetPresenceThreshold=" + probeDatasetPresenceThreshold +
                "\n snpDatasetPresenceThreshold=" + snpDatasetPresenceThreshold +
                "\n probeAndSNPPresenceFilterSampleThreshold=" + probeAndSNPPresenceFilterSampleThreshold +
                "\n runonlypermutation=" + runonlypermutation +
                "\n nrThresds=" + nrThresds +
                "\n probeselection='" + probeselection + '\'' +
                "\n snpselection='" + snpselection + '\'' +
                "\n config=" + config +
                "\n snpprobeselection='" + snpprobeselection + '\'' +
                "\n snpAnnotationFile='" + snpAnnotationFile + '\'' +
                "\n minimalNumberOfDatasets=" + minimalNumberOfDatasets +
                "\n debug=" + debug +
                "\n fullpermutationoutput=" + fullpermutationoutput +
                "\n usetmp=" + usetmp +
                "\n genetosnp='" + genetosnp + '\'' +
                '}';
    }

    private int nrPermutations = 10;
    private int startPermutations = 0;
    private boolean useAbsoluteZscore = false;
    private int finalEQTLBufferMaxLength = 1000000;
    private int nrOfBins = 100;
    private boolean includeSNPsWithoutProperMapping = true;
    private boolean includeProbesWithoutProperMapping = true;

    public Analysis analysisType;

    private int cisdistance = 250000;
    private int transdistance = 5000000;
    private boolean makezscoreplot = true;
    private String probetranslationfile;
    private boolean rescalingOfSampleSize = false;
    private ArrayList<String> featureOccuranceScaleMaps;
    private ArrayList<String> datasetnames;
    private ArrayList<String> datasetPrefix;
    private ArrayList<String> datasetlocations;
    private ArrayList<String> datasetannotations;
    private ArrayList<Integer> selectedProbes;
    private String output;
    private boolean makezscoretable = false;
    private boolean confineSNPs = false;

    private int probeDatasetPresenceThreshold = 0;
    private int snpDatasetPresenceThreshold = 0;
    private int probeAndSNPPresenceFilterSampleThreshold = 0;
    private int runonlypermutation;
    private int nrThresds;
    private String probeselection;
    private String snpselection;
    private XMLConfiguration config;
    private String snpprobeselection;

    private String snpAnnotationFile;
    public int minimalNumberOfDatasets = 1;
    public boolean debug = false;
    public boolean fullpermutationoutput = false;
    public boolean usetmp = false;
    public String genetosnp = null;

    public Analysis getAnalysisType() {
        return analysisType;
    }

    public void copyToOutputDir() throws ConfigurationException {
        config.save(output + "settings.xml");
    }

    public String getGeneToSNP() {
        return genetosnp;
    }

    public void setDatasetPrefix(ArrayList<String> newprefixes) {
        datasetPrefix = newprefixes;
    }


    public enum Analysis {
        CIS, TRANS, CISTRANS
    }

    public void parse(String settings, String settingsTextToReplace, String settingsTextReplaceWith) {
        try {

            if (settingsTextReplaceWith != null) {
                String[] queries = settingsTextToReplace.split(",");
                String[] replacements = settingsTextReplaceWith.split(",");

                if (queries.length != replacements.length) {
                    System.out.println("Issue with replacement strings!");
                    System.out.println("Number of strings to replace not equal!");
                    System.exit(-1);
                }

                for (int q = 0; q < queries.length; q++) {
                    System.out.println("Will replace " + queries[q] + " with " + replacements[q]);
                }

                System.out.println("Will attempt to replace template strings in configuration file.");
                System.out.println(queries.length + " strings to replace with " + replacements.length + " replacements");
                if (replacements.length != queries.length) {
                    System.out.println("Error: number of strings to replace and number of replacements should be equal.");
                    System.exit(-1);
                }

                TextFile tf = new TextFile(settings, TextFile.R);
                String generatedString = RandomStringUtils.randomAlphabetic(12);
//			String generatedString = "v2";//RandomStringUtils.randomAlphabetic(12);
                TextFile tf2 = new TextFile(settings + "-" + generatedString + ".xml", TextFile.W);
                String line = tf.readLine();
                int lnctr = 0;
                while (line != null) {
                    for (int s = 0; s < queries.length; s++) {
                        String query = queries[s];
                        String replacement = replacements[s];
//					System.out.println(lnctr + "\t" + query + "\t" + replacement);
                        if (line.contains(query)) {
//						System.out.println(line + " --> Replacing: " + query + " with " + replacement);
                            line = line.replace(query, replacement);
                        }

                    }
                    tf2.writeln(line);
                    lnctr++;
                    line = tf.readLine();
                }
                tf.close();
                tf2.close();

                config = new XMLConfiguration(settings + "-" + generatedString + ".xml");           // Use the apache XML configuration parser
                Files.delete(Paths.get(settings + "-" + generatedString + ".xml"));
            } else {
                config = new XMLConfiguration(settings);
            }

            nrPermutations = config.getInt("defaults.permutations", 0);
            startPermutations = config.getInt("defaults.startpermutation", 0);
            snpAnnotationFile = config.getString("defaults.snpannotation");
            useAbsoluteZscore = config.getBoolean("defaults.absolutezscore", false);
            finalEQTLBufferMaxLength = config.getInt("defaults.finalnreqtls", 100000);
            cisdistance = config.getInt("defaults.cisprobedistance", 250000);
            transdistance = config.getInt("defaults.transprobedistance", 5000000);
            includeProbesWithoutProperMapping = config.getBoolean("defaults.includeprobeswithoutmapping", true);
            includeSNPsWithoutProperMapping = config.getBoolean("defaults.includesnpswithoutmapping", true);
            makezscoreplot = config.getBoolean("defaults.makezscoreplot", true);
            makezscoretable = config.getBoolean("defaults.makezscoretable", false);
            confineSNPs = config.getBoolean("defaults.confineSNPsToSNPsPresentInAllDatasets", false);
            probetranslationfile = config.getString("defaults.probetranslationfile");

            try {
                usetmp = config.getBoolean("defaults.usetmp", false);
            } catch (Exception e) {

            }
            try {
                fullpermutationoutput = config.getBoolean("defaults.fullpermutationoutput", false);
            } catch (Exception e) {

            }
            try {
                debug = config.getBoolean("defaults.debug", false);
            } catch (Exception e) {

            }

            try {
                minimalNumberOfDatasets = config.getInt("defaults.minimalNumberOfDatasets", 1);
            } catch (Exception e) {

            }

            output = config.getString("defaults.output");
            if (!Gpio.exists(output)) {
                try {
                    Gpio.createDir(output);
                    copyToOutputDir();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }

            probeDatasetPresenceThreshold = config.getInt("defaults.minimalnumberofdatasetsthatcontainprobe", 0);
            snpDatasetPresenceThreshold = config.getInt("defaults.minimalnumberofdatasetsthatcontainsnp", 0);
            probeAndSNPPresenceFilterSampleThreshold = config.getInt("defaults.snpprobeselectsamplesizethreshold", -1);

            runonlypermutation = config.getInt("defaults.runonlypermutation", -1);
            nrThresds = config.getInt("defaults.threads", 1);
            boolean cis = config.getBoolean("defaults.cis", false);
            boolean trans = config.getBoolean("defaults.trans", false);
            if (cis && !trans) {
                analysisType = Analysis.CIS;
            } else if (!cis && trans) {
                analysisType = Analysis.TRANS;
            } else {
                analysisType = Analysis.CISTRANS;
            }

            probeselection = config.getString("defaults.probeselection");

            if (probeselection != null && probeselection.trim().length() == 0) {
                probeselection = null;
            }

            snpselection = config.getString("defaults.snpselection");

            if (snpselection != null && snpselection.trim().length() == 0) {
                snpselection = null;
            }

            snpprobeselection = config.getString("defaults.snpprobeselection");

            if (snpprobeselection != null && snpprobeselection.trim().length() == 0) {
                snpprobeselection = null;
            } else {
                System.out.println("SNP PROBE SELECTION: " + snpprobeselection);
            }

//			try {
            genetosnp = config.getString("defaults.genetosnp");

//			} catch (NoSuchElementException e) {
//
//			}

            int i = 0;

            String dataset = "";
            datasetnames = new ArrayList<String>();
            datasetlocations = new ArrayList<String>();
            datasetannotations = new ArrayList<String>();
            datasetPrefix = new ArrayList<String>();
            featureOccuranceScaleMaps = new ArrayList<String>();

            while (dataset != null) {
                dataset = config.getString("datasets.dataset(" + i + ").name");  // see if a dataset is defined
                if (dataset != null) {

                    datasetnames.add(dataset);
                    String prefix = config.getString("datasets.dataset(" + i + ").prefix");  // see if a dataset is defined

                    if (prefix == null) {
                        prefix = "Dataset";
                    }
                    datasetPrefix.add(prefix);
                    String datasetlocation = config.getString("datasets.dataset(" + i + ").location");  // see if a dataset is defined
                    String datasetannotation = config.getString("datasets.dataset(" + i + ").expressionplatform");  // see if a dataset is defined

                    datasetlocations.add(datasetlocation);
                    datasetannotations.add(datasetannotation);
                    String featureOccuranceScaleMap = config.getString("datasets.dataset(" + i + ").featureOccuranceScaleMap");
                    if (featureOccuranceScaleMap != null && featureOccuranceScaleMap.trim().length() == 0) {
                        featureOccuranceScaleMap = null;
                    }
                    if (featureOccuranceScaleMap != null) {
                        System.out.println("Feature rescaling values: " + featureOccuranceScaleMap);
                        rescalingOfSampleSize = true;
                    }
                    featureOccuranceScaleMaps.add(featureOccuranceScaleMap);
                }
                i++;
            }

            summarize();
            // parse datasets
        } catch (ConfigurationException | IOException e) {
            e.printStackTrace();
        }
    }

    private void summarize() {
        System.out.println(this.toString());
    }

    public int getStartPermutations() {
        return startPermutations;
    }

    public void setStartPermutations(int startPermutations) {
        this.startPermutations = startPermutations;
    }

    public boolean isConfineSNPs() {
        return confineSNPs;
    }

    public void setConfineSNPs(boolean confineSNPs) {
        this.confineSNPs = confineSNPs;
    }

    /**
     * @return the nrPermutations
     */
    public int getNrPermutations() {
        return nrPermutations;
    }

    /**
     * @param nrPermutations the nrPermutations to set
     */
    public void setNrPermutations(int nrPermutations) {
        this.nrPermutations = nrPermutations;
    }

    /**
     * @return the useAbsoluteZscore
     */
    public boolean isUseAbsoluteZscore() {
        return useAbsoluteZscore;
    }

    /**
     * @param useAbsoluteZscore the useAbsoluteZscore to set
     */
    public void setUseAbsoluteZscore(boolean useAbsoluteZscore) {
        this.useAbsoluteZscore = useAbsoluteZscore;
    }

    /**
     * @return the finalEQTLBufferMaxLength
     */
    public int getFinalEQTLBufferMaxLength() {
        return finalEQTLBufferMaxLength;
    }

    /**
     * @param finalEQTLBufferMaxLength the finalEQTLBufferMaxLength to set
     */
    public void setFinalEQTLBufferMaxLength(int finalEQTLBufferMaxLength) {
        this.finalEQTLBufferMaxLength = finalEQTLBufferMaxLength;
    }

    /**
     * @return the nrOfBins
     */
    public int getNrOfBins() {
        return nrOfBins;
    }

    /**
     * @param nrOfBins the nrOfBins to set
     */
    public void setNrOfBins(int nrOfBins) {
        this.nrOfBins = nrOfBins;
    }

    /**
     * @return the includeSNPsWithoutProperMapping
     */
    public boolean isIncludeSNPsWithoutProperMapping() {
        return includeSNPsWithoutProperMapping;
    }

    /**
     * @param includeSNPsWithoutProperMapping the
     *                                        includeSNPsWithoutProperMapping to set
     */
    public void setIncludeSNPsWithoutProperMapping(boolean includeSNPsWithoutProperMapping) {
        this.includeSNPsWithoutProperMapping = includeSNPsWithoutProperMapping;
    }

    /**
     * @return the includeProbesWithoutProperMapping
     */
    public boolean isIncludeProbesWithoutProperMapping() {
        return includeProbesWithoutProperMapping;
    }

    /**
     * @param includeProbesWithoutProperMapping the
     *                                          includeProbesWithoutProperMapping to set
     */
    public void setIncludeProbesWithoutProperMapping(boolean includeProbesWithoutProperMapping) {
        this.includeProbesWithoutProperMapping = includeProbesWithoutProperMapping;
    }

    /**
     * @return the cisdistance
     */
    public int getCisdistance() {
        return cisdistance;
    }

    /**
     * @param cisdistance the cisdistance to set
     */
    public void setCisdistance(int cisdistance) {
        this.cisdistance = cisdistance;
    }

    /**
     * @return the transdistance
     */
    public int getTransdistance() {
        return transdistance;
    }

    /**
     * @param transdistance the transdistance to set
     */
    public void setTransdistance(int transdistance) {
        this.transdistance = transdistance;
    }

    /**
     * @return the makezscoreplot
     */
    public boolean isMakezscoreplot() {
        return makezscoreplot;
    }

    /**
     * @param makezscoreplot the makezscoreplot to set
     */
    public void setMakezscoreplot(boolean makezscoreplot) {
        this.makezscoreplot = makezscoreplot;
    }

    /**
     * @return the probetranslationfile
     */
    public String getProbetranslationfile() {
        return probetranslationfile;
    }

    /**
     * @param probetranslationfile the probetranslationfile to set
     */
    public void setProbetranslationfile(String probetranslationfile) {
        this.probetranslationfile = probetranslationfile;
    }

    /**
     * @return the datasetnames
     */
    public ArrayList<String> getDatasetnames() {
        return datasetnames;
    }

    /**
     * @param datasetnames the datasetnames to set
     */
    public void setDatasetnames(ArrayList<String> datasetnames) {
        this.datasetnames = datasetnames;
    }

    /**
     * @return the datasetlocations
     */
    public ArrayList<String> getDatasetlocations() {
        return datasetlocations;
    }

    /**
     * @param datasetlocations the datasetlocations to set
     */
    public void setDatasetlocations(ArrayList<String> datasetlocations) {
        this.datasetlocations = datasetlocations;
    }

    /**
     * @return the datasetannotations
     */
    public ArrayList<String> getDatasetannotations() {
        return datasetannotations;
    }

    /**
     * @param featureOccuranceScaleMaps the featureOccuranceScaleMaps to set
     */
    public void setFeatureOccuranceScaleMaps(ArrayList<String> featureOccuranceScaleMaps) {
        this.featureOccuranceScaleMaps = featureOccuranceScaleMaps;
    }

    /**
     * @return the featureOccuranceScaleMaps
     */
    public ArrayList<String> getFeatureOccuranceScaleMaps() {
        return featureOccuranceScaleMaps;
    }

    /**
     * @return the rescalingOfSampleSize
     */
    public boolean getRescalingOfSampleSize() {
        return rescalingOfSampleSize;
    }

    /**
     * @param datasetannotations the datasetannotations to set
     */
    public void setDatasetannotations(ArrayList<String> datasetannotations) {
        this.datasetannotations = datasetannotations;
    }

    /**
     * @return the output
     */
    public String getOutput() {
        return output;
    }

    /**
     * @param output the output to set
     */
    public void setOutput(String output) {
        this.output = output;
    }

    /**
     * @return the makezscoretable
     */
    public boolean isMakezscoretable() {
        return makezscoretable;
    }

    /**
     * @param makezscoretable the makezscoretable to set
     */
    public void setMakezscoretable(boolean makezscoretable) {
        this.makezscoretable = makezscoretable;
    }

    /**
     * @return the probeDatasetPresenceThreshold
     */
    public int getProbeDatasetPresenceThreshold() {
        return probeDatasetPresenceThreshold;
    }

    /**
     * @param probeDatasetPresenceThreshold the probeDatasetPresenceThreshold to
     *                                      set
     */
    public void setProbeDatasetPresenceThreshold(int probeDatasetPresenceThreshold) {
        this.probeDatasetPresenceThreshold = probeDatasetPresenceThreshold;
    }

    /**
     * @return the snpDatasetPresenceThreshold
     */
    public int getSnpDatasetPresenceThreshold() {
        return snpDatasetPresenceThreshold;
    }

    /**
     * @param snpDatasetPresenceThreshold the snpDatasetPresenceThreshold to set
     */
    public void setSnpDatasetPresenceThreshold(int snpDatasetPresenceThreshold) {
        this.snpDatasetPresenceThreshold = snpDatasetPresenceThreshold;
    }

    /**
     * @return the probeAndSNPPresenceFilterSampleThreshold
     */
    public int getProbeAndSNPPresenceFilterSampleThreshold() {
        return probeAndSNPPresenceFilterSampleThreshold;
    }

    /**
     * @param probeAndSNPPresenceFilterSampleThreshold the
     *                                                 probeAndSNPPresenceFilterSampleThreshold to set
     */
    public void setProbeAndSNPPresenceFilterSampleThreshold(int probeAndSNPPresenceFilterSampleThreshold) {
        this.probeAndSNPPresenceFilterSampleThreshold = probeAndSNPPresenceFilterSampleThreshold;
    }

    /**
     * @return the runonlypermutation
     */
    public int getRunonlypermutation() {
        return runonlypermutation;
    }

    /**
     * @param runonlypermutation the runonlypermutation to set
     */
    public void setRunonlypermutation(int runonlypermutation) {
        this.runonlypermutation = runonlypermutation;
    }

    /**
     * @return the nrThresds
     */
    public int getNrThreads() {
        return nrThresds;
    }

    /**
     * @param nrThresds the nrThresds to set
     */
    public void setNrThresds(int nrThresds) {
        this.nrThresds = nrThresds;
    }

    public ArrayList<String> getDatasetPrefix() {
        return datasetPrefix;
    }

    /**
     * @return the probeselection
     */
    public String getProbeselection() {
        return probeselection;
    }

    /**
     * @param probeselection the probeselection to set
     */
    public void setProbeselection(String probeselection) {
        this.probeselection = probeselection;
    }

    public String getSNPSelection() {
        return snpselection;
    }

    public String getSNPProbeSelection() {
        return snpprobeselection;
    }

    void save() {
        try {
            config.save(output + "metasettings.xml");
        } catch (ConfigurationException ex) {
            throw new RuntimeException(ex);
        }

    }

    public String getSNPAnnotationFile() {
        return snpAnnotationFile;
    }


    public XMLConfiguration getConfig() {
        return config;
    }

    public void save(String loc) throws ConfigurationException {
        config.save(loc);
    }
}
