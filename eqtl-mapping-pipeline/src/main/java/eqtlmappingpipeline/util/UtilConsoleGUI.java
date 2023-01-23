/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.util;

import eqtlmappingpipeline.binarymeta.Main;
import eqtlmappingpipeline.metaqtl3.FDR;
import eqtlmappingpipeline.metaqtl3.FDR.FDRMethod;
import eqtlmappingpipeline.pcaoptimum.PCAOptimum;
import eqtlmappingpipeline.textmeta.FixedEffectMetaAnalysis;
import umcg.genetica.console.ConsoleGUIElems;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.trityper.ConvertDoubleMatrixDataToTriTyper;
import umcg.genetica.io.trityper.util.ChrAnnotation;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;


import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.UnsupportedEncodingException;
import java.util.Arrays;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * @author harmjan
 */
public class UtilConsoleGUI {

    public static enum MODE {

        GETSNPSFROMREGION, GETSNPSINPROBEREGION, FDR, GETMAF, MERGE, REGRESS, GETSNPSTATS, PROXYSEARCH, DOTPLOT, META,
        SORTEQTLFILEBYZSCORE, CONVERTTRITYPERTOMATRIX, CONVERTBINARYMATRIX, GETSNPPROBECOMBINATIONS, NONGENETICPCACORRECTION, REGRESSKNOWN, CREATTTFROMDOUBLEMAT,
        ADDANNOTATIONTOQTLFILE, LOOKUPEFFECTS, FDRPROBE, PHENOTYPESAMPLEFILTER, SPLITEQTLFILEBYCHR, QTLFILEMERGE, EQTLEQTMLINK, SPLITPHENO,
        SORTEQTLFILEBYPVALUE, SORTEQTLFILEBYCHRPOS, SPLITTRITYPERBYCHR, GETMAFFROMQCLOG, CALCULATEBETA, CONVERTTOSMR, CONVERTBINMAT, QQPLOT
    }

    ;
    MODE run;

    public UtilConsoleGUI(String[] args) {

        String settingstexttoreplace = null;
        String settingstexttoreplacewith = null;
//        boolean cis = false;
//        boolean trans = false;
//        String outtype = "text";
//        String inexpplatform = null;
//        Integer threads = null;
        String settingsfile = null;
        String in = null;
        String in2 = null;
        String out = null;
        int perm = 0;
        String inexp = null;
        String inexpannot = null;
        String gte = null;
        String snpfile = null;
        String probefile = null;
        String region = "";
        String ranked = "";

        String annot = null;
        String probeselectionlist = null;
        Integer stepSize = 5;
        Integer max = 5;
        String fileQtlsToRegressOut = null;

        Double threshold = null;
        Integer nreqtls = null;

        Double r2 = null;
        Double maf = 0.05;
        Double cr = 0.95;
        Double hwep = 0.001;
        Integer dist = 1000000;
        Integer threads = 1;
        Integer minnrdatasets = null;
        Integer minnrsamples = null;

        String snpprobeselectionlist = null;
        boolean createQQPlot = true;
        boolean createLargeFdrFile = true;
        boolean stringentFDR = false;
        boolean sortsnps = true;
        boolean vcf = false;
        boolean sortById = false;
        boolean trityper = false;
        String sources = null;
        String keyValuePairs = null;
        String annotationIds = null;
        String geneAnnotationFile = null;

        FDRMethod fdrMethod = FDRMethod.ALL;

        for (int i = 0; i < args.length; i++) {
            String arg = args[i];
            String val = null;

            if (i + 1 < args.length) {
                val = args[i + 1];
            }

            if (arg.equals("--tritypertomatrix")) {
                region = val;
                run = MODE.CONVERTTRITYPERTOMATRIX;
            } else if (arg.equals("--convertbinarymatrix")) {
                region = val;
                run = MODE.CONVERTBINARYMATRIX;
            } else if (arg.equals("--eqtmlink")) {
                run = MODE.EQTLEQTMLINK;
            } else if (arg.equals("--mergeqtlfile")) {
                run = MODE.QTLFILEMERGE;
            } else if (arg.equals("--splitpheno")) {
                run = MODE.SPLITPHENO;
            } else if (arg.equals("--spliteqtlfilebychr")) {
                run = MODE.SPLITEQTLFILEBYCHR;
            } else if (arg.equals("--splitTT")) {
                run = MODE.SPLITTRITYPERBYCHR;
            } else if (arg.equals("--filterpheno")) {
                run = MODE.PHENOTYPESAMPLEFILTER;
            } else if (arg.equals("--getsnpsinregion")) {
                region = val;
                run = MODE.GETSNPSFROMREGION;
            } else if (arg.equals("--sorteqtlfilebyz") || arg.equals("--sortfile")) {
                region = val;
                run = MODE.SORTEQTLFILEBYZSCORE;
            } else if (arg.equals("--sorteqtlfilebyp")) {
                region = val;
                run = MODE.SORTEQTLFILEBYPVALUE;
            } else if (arg.equals("--sorteqtlfilebypos")) {
                region = val;
                run = MODE.SORTEQTLFILEBYCHRPOS;
            } else if (arg.equals("--calculatebeta")) {
                run = MODE.CALCULATEBETA;
            } else if (arg.equals("--findproxy")) {
                region = val;
                run = MODE.PROXYSEARCH;
            } else if (arg.equals("--getmaf")) {
                region = val;
                run = MODE.GETMAF;
            } else if (arg.equals("--getmaffromqclog")) {
                run = MODE.GETMAFFROMQCLOG;
            } else if (arg.equals("--getsnpsinproberegion")) {
                region = val;
                run = MODE.GETSNPSINPROBEREGION;
            } else if (arg.equals("--merge")) {
                run = MODE.MERGE;
            } else if (arg.equals("--fdr")) {
                region = val;
                run = MODE.FDR;
            } else if (arg.equals("--fdrprobe")) {
                region = val;
                run = MODE.FDRPROBE;
            } else if (arg.equals("--dotplot")) {
                region = val;
                run = MODE.DOTPLOT;
            } else if (arg.equals("--regress")) {
                run = MODE.REGRESS;
            } else if (arg.equals("--snpstats")) {
                run = MODE.GETSNPSTATS;
            } else if (arg.equals("--meta")) {
                run = MODE.META;
            } else if (arg.equals("--regressknown")) {
                run = MODE.REGRESSKNOWN;
            } else if (arg.equals("--getSNPProbeCombinations")) {
                run = MODE.GETSNPPROBECOMBINATIONS;
            } else if (arg.equals("--nonGeneticPcaCorrection")) {
                run = MODE.NONGENETICPCACORRECTION;
            } else if (arg.equals("--formatAsTT")) {
                run = MODE.CREATTTFROMDOUBLEMAT;
            } else if (arg.equals("--settings")) {
                settingsfile = val;
            } else if (arg.equals("--in")) {
                in = val;
            } else if (arg.equals("--sources")) {
                sources = val;
            } else if (arg.equals("--keyValuePairs")) {
                keyValuePairs = val;
            } else if (arg.equals("--idsToAnnotate")) {
                annotationIds = val;
            } else if (arg.equals("--geneAnnotation")) {
                geneAnnotationFile = val;
            } else if (arg.equals("--in2")) {
                in2 = val;
            } else if (arg.equals("--out")) {
                out = val;
            } else if (arg.equals("--inexp")) {
                inexp = val;
            } else if (arg.equals("--inexpannot")) {
                inexpannot = val;
            } else if (arg.equals("--gte")) {
                gte = val;
            } else if (args[i].equals("--annot")) {
                annot = val;
            } else if (args[i].equals("--fdrmethod")) {
                val = val.toLowerCase();
                if (val.equals("probe")) {
                    fdrMethod = FDRMethod.PROBELEVEL;
                } else if (val.equals("gene")) {
                    fdrMethod = FDRMethod.GENELEVEL;
                } else if (val.equals("snp")) {
                    fdrMethod = FDRMethod.SNPLEVEL;
                } else if (val.equals("full")) {
                    fdrMethod = FDRMethod.FULL;
                }
            } else if (args[i].equals("--stringentFDR")) {
                stringentFDR = true;
            } else if (arg.equals("--snps")) {
                snpfile = val;
            } else if (arg.equals("--probes")) {
                probefile = val;
            } else if (arg.equals("--threads")) {
                threads = Integer.parseInt(val);
            } else if (arg.equals("--perm")) {
                perm = Integer.parseInt(val);
            } else if (arg.equals("--nreqtls")) {
                nreqtls = Integer.parseInt(val);
            } else if (arg.equals("--dontsortsnps")) {
                sortsnps = false;
            } else if (arg.equals("--threshold")) {
                threshold = Double.parseDouble(val);
            } else if (arg.equals("--r2")) {
                r2 = Double.parseDouble(val);
            } else if (arg.equals("--maf")) {
                maf = Double.parseDouble(val);
            } else if (arg.equals("--hwep")) {
                hwep = Double.parseDouble(val);
            } else if (arg.equals("--dist")) {
                dist = Integer.parseInt(val);
            } else if (arg.equals("--skipqqplot")) {
                createQQPlot = false;
            } else if (arg.equals("--skipLargeFDRFile")) {
                createLargeFdrFile = false;
            } else if (args[i].equals("--snpselectionlist")) {
                snpfile = val;
            } else if (args[i].equals("--vcf")) {
                vcf = true;
            } else if (args[i].equals("--sortbyid")) {
                sortById = true;
            } else if (args[i].equals("--trityper")) {
                trityper = true;
            } else if (args[i].equals("--probeselectionlist")) {
                probeselectionlist = val;
            } else if (args[i].equals("--snpprobeselectionlist")) {
                snpprobeselectionlist = val;
            } else if (args[i].equals("--stepsizepcaremoval")) {
                stepSize = Integer.parseInt(val);
            } else if (args[i].equals("--maxnrpcaremoved")) {
                max = Integer.parseInt(val);
            } else if (args[i].equals("--QTLS")) {
                fileQtlsToRegressOut = val;
            } else if (args[i].equals("--preRankExpressionFiles")) {
                ranked = "ranked";
            } else if (arg.equals("--replacetext")) {
                settingstexttoreplace = val;
            } else if (arg.equals("--replacetextwith")) {
                settingstexttoreplacewith = val;
            } else if (arg.equals("--converteqtlfiletosmr")) {
                run = MODE.CONVERTTOSMR;
            } else if (arg.equals("--qqplot")) {
                run = MODE.QQPLOT;
            }
//            }  else if (arg.equals("--inexpplatform")) {
//                inexpplatform = val;
//            }

        }
        if (run == null) {
            System.err.println("Please specify an util.");
            printUsage();
        } else {
            try {
                switch (run) {
                    case CONVERTTRITYPERTOMATRIX:
                        if (settingsfile == null) {
                            System.out.println("Usage: --tritypertomatrix --settings settings.xml [--vcf|--trityper] [--sortbyid] ");
                        } else {
                            TriTyperToDosageMatrix ttd = new TriTyperToDosageMatrix();
                            if (trityper) {
                                ttd.runTriTyper(settingsfile, settingstexttoreplace, settingstexttoreplacewith, sortById);
                            } else {
                                ttd.run(settingsfile, settingstexttoreplace, settingstexttoreplacewith, vcf, sortById);
                            }
                        }
                        break;
                    case CONVERTBINARYMATRIX:
                        if (in == null || out == null) {
                            System.out.println("Usage: --util --convertbinarymatrix --in /path/to/matrix.binary --out /path/to/textoutput.txt");
                        } else {
                            if (in.endsWith(".txt")) {
                                System.out.println("The file provided with --in is already a text file: " + in);
                            } else {
                                if (in.endsWith(".dat")) {
                                    in = in.substring(0, in.length() - 4);
                                }
                                System.out.println("Converting: " + in);
                                DoubleMatrixDataset<String, String> ds = DoubleMatrixDataset.loadDoubleBinaryData(in);
                                ds.save(out);
                            }
                        }
                        break;
                    case SPLITPHENO:
                        BatchSplitter bs = new BatchSplitter();
                        if (in == null) {
                            System.out.println("Use --in");
                        }
                        if (inexpannot == null) {
                            System.out.println("Use --inexpannot");
                        }
                        if (out == null) {
                            System.out.println("Use --out");
                        }
                        bs.splitPhenotypePerChr(in, inexpannot, out);
                        break;
                    case SPLITTRITYPERBYCHR:
                        TriTyperSplitter t = new TriTyperSplitter();
                        t.split(in, out);
                        break;
                    case REGRESS:

                        RegressCisEffectsFromGeneExpressionData r = new RegressCisEffectsFromGeneExpressionData(args);
                        break;

                    case PROXYSEARCH:

                        if (in == null || snpfile == null || out == null || r2 == null) {
                            System.out.println("Usage: --mode util --findproxy --r2 0.8 --snps snpfile.txt --out outfile --in /Path/To/TriTyperReference/ [--hwep 0.001] [--maf 0.05] [--cr 0.95]");
                        } else {
                            LDCalculator.proxyLookUpInReferenceDataset(in, snpfile, maf, hwep, cr, r2, out, dist);
                        }
                        break;

                    case MERGE:

                        if (in == null || region == null) {
                            System.out.println("USAGE: --merge --in dataset --in2 dataset2 --out outdir [--snps snpfile]");
                        } else {
                            GenotypeDataMerger m = new GenotypeDataMerger();
                            m.merge(in, in2, out, snpfile);
                        }
                        break;

                    case GETMAF:

                        if (in == null || region == null) {
                            System.out.println("USAGE: --getmaf snplistfile --in dataset");
                        } else {
                            GenotypeDataQuery dq = new GenotypeDataQuery();
                            dq.getSNPMAF(in, region);
                        }
                        break;
                    case GETMAFFROMQCLOG:
                        if (in == null || out == null) {
                            System.out.println("USAGE: --getmaffromqclog --in snpqclog --out out");
                        } else {
                            GetMAFFromSNPQCLog q = new GetMAFFromSNPQCLog();
                            q.run(in, out);
                        }
                        break;

                    case GETSNPSTATS:

                        if (in == null || region == null) {
                            System.out.println("USAGE: --in dataset");
                        } else {
                            GenotypeDataQuery dq = new GenotypeDataQuery();
                            if (in2 != null) {
                                dq.getSNPStatsForAllSNPs(in, in2);
                            } else {
                                dq.getSNPStatsForAllSNPs(in);
                            }
                        }
                        break;
                    case QTLFILEMERGE:
                        if (in == null || out == null) {
                            System.out.println("USAGE: --in eQTLFile --out eQTLFile\nSpecify multiple files using CHR template, or provide comma separated list.");
                        } else {
                            QTLFileMerger m = new QTLFileMerger();
                            m.mergeChr(in, out);
                        }
                        break;
                    case SORTEQTLFILEBYPVALUE:
                        if (in == null) {
                            System.out.println("USAGE: --in eQTLFile --out eQTLFile");
                        } else {
                            QTLFileSorter f = new QTLFileSorter();
                            int n = 2500000;
                            if (nreqtls != null) {
                                n = nreqtls;
                            }
                            f.run(in, out, n, QTLFileSorter.SORTBY.P);
                        }
                        break;
                    case SORTEQTLFILEBYZSCORE:
                        if (in == null) {
                            System.out.println("USAGE: --in eQTLFile --out eQTLFile");
                        } else {
                            QTLFileSorter f = new QTLFileSorter();
                            int n = 2500000;
                            if (nreqtls != null) {
                                n = nreqtls;
                            }
                            f.run(in, out, n, QTLFileSorter.SORTBY.Z);
                        }
                        break;
                    case SORTEQTLFILEBYCHRPOS:
                        if (in == null) {
                            System.out.println("USAGE: --in eQTLFile --out eQTLFile");
                        } else {
                            QTLFileSorter f = new QTLFileSorter();
                            int n = 2500000;
                            if (nreqtls != null) {
                                n = nreqtls;
                            }
                            f.run(in, out, n, QTLFileSorter.SORTBY.POS);
                        }
                        break;
                    case SPLITEQTLFILEBYCHR:
                        if (in == null) {
                            System.out.println("USAGE: --in eQTLFile --out eQTLFileprefix [--nreqtls 2500000]");
                        } else {
                            QTLFileSplitter s = new QTLFileSplitter();
                            s.split(in, out);
                        }
                        break;
                    case GETSNPSFROMREGION:
                        if (in == null || region == null) {
                            System.out.println("To use --getsnpsfromregion, please use --in to point to the genotype data and supply a region to query.");
                            printUsage();
                        } else {
                            int chr = -1;
                            int chrposA = -1;
                            int chrposB = -1;
                            GenotypeDataQuery q = new GenotypeDataQuery();
                            try {
                                String[] elems = region.split(":");
                                chr = ChrAnnotation.parseChr(elems[0]);
                                elems = elems[1].split("-");
                                chrposA = Integer.parseInt(elems[0]);
                                chrposB = Integer.parseInt(elems[1]);
                            } catch (Exception e) {
                                System.err.println("Error: malformed query: " + region);
                                ;
                            }
                            q.getSNPsInRegion(in, chr, chrposA, chrposB);
                        }

                        break;

                    case GETSNPSINPROBEREGION:
                        if (snpfile == null || inexpannot == null || probefile == null) {
                            System.out.println("To use --getsnpsinproberegion, please use --snps, --probes, and --inexpannot");
                            printUsage();
                        } else {
                            ProbeSNPMapper psm = new ProbeSNPMapper();
                            psm.mapprobes(snpfile, inexpannot, probefile);
                        }

                        break;

                    case FDR:
                        if (in == null || threshold == null || nreqtls == null || perm == 0) {
                            System.out.println("To use --fdr, please use --in, --threshold, and --perm and --nreqtls");
                            System.out.println("Optional: --snpselectionlist, --probeselectionlist, --snpprobeselectionlist");
                            printUsage();
                        } else {
                            if (snpfile != null || snpprobeselectionlist != null || probeselectionlist != null) {
                                try {
                                    FDR.calculateFDRAdvanced(in, perm, nreqtls, threshold, createQQPlot, null, null, fdrMethod, createLargeFdrFile, snpfile, probeselectionlist, snpprobeselectionlist);
                                } catch (IOException e) {
                                    e.printStackTrace();
                                    System.exit(1);
                                }
                            } else {
                                try {
                                    FDR.calculateFDR(in, perm, nreqtls, threshold, createQQPlot, null, null, fdrMethod, createLargeFdrFile);
                                } catch (IOException e) {
                                    e.printStackTrace();
                                    System.exit(1);
                                }
                            }
                        }

                        break;
                    case FDRPROBE:
                        if (in == null || threshold == null || perm == 0) {
                            System.out.println("To use --fdrprobe, please use --in, --threshold, and --perm");
                            System.out.println("Optional: --stringentFDR");
                            printUsage();
                        } else {
                            try {
                                ProbeSpecificFDR.calculateFDR(in, perm, threshold, false, stringentFDR, createLargeFdrFile);
                            } catch (IOException e) {
                                e.printStackTrace();
                                System.exit(1);
                            }
                        }

                        break;

                    case META:
                        if (in == null || out == null) {
                            System.out.println("To use --meta, please use --in, and --out");
                            printUsage();
                        } else {
                            FixedEffectMetaAnalysis f = new FixedEffectMetaAnalysis();
                            f.run(in, out, minnrdatasets, minnrsamples);
                        }

                        break;

                    case DOTPLOT:
                        if (in == null) {
                            System.out.println("Usage: --dotplot --in /path/to/file.txt");
                        } else {
                            QTLDotPlotter d = new QTLDotPlotter();
                            d.plot(in);
                        }
                        break;

                    case GETSNPPROBECOMBINATIONS:

                        try {
                            NoLdSnpProbeListCreator.main(Arrays.copyOfRange(args, 2, args.length));
                        } catch (UnsupportedEncodingException ex) {
                            throw new RuntimeException(ex);
                        } catch (FileNotFoundException ex) {
                            throw new RuntimeException(ex);
                        } catch (IOException ex) {
                            throw new RuntimeException(ex);
                        } catch (Exception ex) {
                            throw new RuntimeException(ex);
                        }

                        break;

                    case NONGENETICPCACORRECTION:

                        boolean cmdlinesettings = (in == null || out == null || inexp == null || gte == null);

                        if (settingsfile == null && cmdlinesettings) {
                            System.out.println("Please specify --settings, --in, --out, --stepsizepcaremoval, --maxnrpcaremoved, --gte, --nreqtls and --dontsortsnps");
                        } else {
                            try {


                                PCAOptimum p = new PCAOptimum();
//            public void alternativeInitialize(String ingt, String inexp, String inexpplatform, String inexpannot, String gte, String out, boolean cis, boolean trans, int perm, String snpfile, Integer threads) throws IOException, Exception {


                                perm = 10;
                                p.alternativeInitialize(settingsfile, in, inexp, null, annot, gte, out, true, true, perm, snpfile, threads, sortsnps);
                                p.performeQTLMappingOverEigenvectorMatrixAndReNormalize(out, stepSize, max, nreqtls);
                            } catch (IOException ex) {
                                throw new RuntimeException(ex);
                            } catch (Exception ex) {
                                throw new RuntimeException(ex);
                            }

                        }
                        break;

                    case REGRESSKNOWN:
                        if (settingsfile == null || fileQtlsToRegressOut == null) {
                            System.out.println("Please specify --settings, --QTLS");
                            break;
                        }

                        if (!Gpio.exists(fileQtlsToRegressOut)) {
                            System.err.println("ERROR: you have specified an eQTL file to regress out, but the file was not found " + fileQtlsToRegressOut);
                            break;
                        }
                        RegressCisEffectsFromGeneExpressionData regress = new RegressCisEffectsFromGeneExpressionData(settingsfile, fileQtlsToRegressOut, settingstexttoreplace, settingstexttoreplacewith);
                        break;

                    case CREATTTFROMDOUBLEMAT:
                        if (inexpannot == null || in == null || out == null) {
                            System.out.println("Please specify --inexpannot, --in, --out");
                        } else {
                            String[] argsNew = {"-m", inexpannot, "-d", in, "-o", out, "-r", ranked};

                            ConvertDoubleMatrixDataToTriTyper.main(argsNew);
                        }
                        break;

                    case ADDANNOTATIONTOQTLFILE:
                        QTLAnnotator.addAnnotationToQTLOutput(in, sources, keyValuePairs, annotationIds, geneAnnotationFile, out);
                        break;

                    case LOOKUPEFFECTS:
                        if (in2 == null || in == null || out == null) {
                            System.out.println("Please specify --in, --in2, --out");
                        } else {
                            QTLLookup.lookUpEffects(in, in2, out);
                        }
                        break;
                    case PHENOTYPESAMPLEFILTER:
                        if (in == null || out == null || annot == null) {
                            System.out.println("Please specify --in, --out, --annot");
                        } else {
                            ExpressionSampleFilter f = new ExpressionSampleFilter();
                            f.filter(in, out, annot);
                        }
                        break;
                    case CALCULATEBETA:
                        if (in == null || out == null || in2 == null) {
                            System.out.println("Please specify --in eqtlfile, --out outeqtlfile, --in2 maffile");
                        } else {
                            ReplaceBeta rpb = new ReplaceBeta();

                            rpb.run(in, out, in2);
                        }
                        break;
                    case CONVERTTOSMR:
                        if (in == null || in2 == null || out == null) {
                            System.out.println("Usage: --in eqtlfile --out eqtlfile --in2 snpqcfile");
                        } else {
                            ConvertTOSMRFormat c = new ConvertTOSMRFormat();

                            c.run(in, in2, out);
                        }
                        break;
                    case QQPLOT:
                        if (in == null || threshold == null || nreqtls == null || perm == 0) {
                            System.out.println("To use --fdr, please use --in, --threshold, and --perm and --nreqtls");
                            System.out.println("Optional: --snpselectionlist, --probeselectionlist, --snpprobeselectionlist");
                            printUsage();
                        } else {
                            QQPlotter q = new QQPlotter();
                            q.run(in, out, threshold, perm, nreqtls);
                        }
                        break;

                }
            } catch (Exception e) {
                e.printStackTrace();
                System.exit(-1);
            }
        }
    }

    private void printUsage() {
        System.out.print("\tUtil\n" + ConsoleGUIElems.LINE);
        System.out.println("Util contains small utilities.");

        System.out.println("");
        System.out.print("Available Utilities:\n" + ConsoleGUIElems.LINE);

        System.out.println("--getsnpsinregion\t\tGet SNPs in a certain region: chr positionA positionB: Y:12000-13000 would get all SNPs on chr Y between 12000 and 13000 bp\n"
                + "--getsnpsinproberegion\t\tGet SNPs in a certain set of probes (specify with --probes)\n"
                + "--fdr\t\t\t\tCalculated FDR.\n"
                + "--fdrprobe\t\t\t\tCalculated FDR per probe.\n"
                + "--getmaf\t\t\tGets maf for snp\n"
                + "--merge\t\t\t\tMerges two datasets\n"
                + "--snpstats\t\t\tGets HWE, MAF, and CR for all SNPs\n"
                + "--findproxy\t\t\tSearches for a proxy given a list of SNPs\n"
                + "--dotplot\t\t\tCreates dotplot from eQTL result file\n"
                + "--regress\t\t\tRemoves eQTL effects from gene expression data.\n"
                + "--regressknown\t\t\tRemoves known cis-eQTL effects from gene expression data.\n"
                + "--sortfile\t\t\t(or --sorteqtlfilebyz) Sort eQTL files by abs(Z-score).\n"
                + "--sorteqtlfilebypos\t\t\tSort eQTL files by chromosome position.\n"
                + "--calculatebeta\t\t\tCalculate beta using Z-score, MAF, and N\n"
                + "--getmaffromqclog\t\t\tCalculate MAF from SNPQCLog.txt.gz, for use with --calculatebeta\n"
                + "--meta\t\t\t\tFixed effect meta analysis.\n"
                + "--nonGeneticPcaCorrection\tCorrect expression data for non-genetic components.\n"
                + "--getSNPProbeCombinatios\tCreate list of valid SNP-Probe combinations to test.\n"
                + "--formatAsTT\t\t\tConverts a doublematrix dataset to a TriTyper genotype file.\n"
                + "--convertbinarymatrix\t\tConverts binary matrix to text\n"
                + "--filterpheno\t\tFilters phenotype file for sample list\n"
                + "--splitpheno\t\tSplits phenotype file by chromosome\n"
                + "--splitTT\t\tSplit trityper folder into chromosome chunks\n"
                + "--spliteqtlfilebychr\t\tSplit eQTL file by chromosome\n"
                + "--mergeqtlfile\t\tMerge QTL files (and sort them)\n"
                + "--eqtmlink\t\tLink eQTM and eQTL files based on probe/gene name\n"
                + "--fdrmethod\t\tEither probe, gene, snp or full\n"
                + "--converteqtlfiletosmr\t\tConvert eQTL file to SMR format\n"
                + "--tritypertomatrix\t\tConvert (set of) TriTyper files to a text-based dosage matrix. Use --vcf or --trityper to output in vcf or TriTyper formats\n"
                + "--qqplot\t\tCreate QQ plot from eQTL output dir");
        System.out.println("");

    }
}
