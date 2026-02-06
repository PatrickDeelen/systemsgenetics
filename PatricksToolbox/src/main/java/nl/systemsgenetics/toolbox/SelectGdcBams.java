package nl.systemsgenetics.toolbox;

import com.opencsv.CSVWriter;
import org.json.JSONArray;
import org.json.JSONObject;

import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;

public class SelectGdcBams {


    public static void main(String[] args) throws IOException {

        // Read entire file into a String
        String jsonContent = Files.readString(Path.of("C:\\Users\\patri\\OneDrive - UMCG\\Documents - Functional Genomics\\Projects\\GutEqtls\\GDC\\metadata.repository.2026-01-26-dna-bam.json"));

        JSONArray jsonArray = new JSONArray(jsonContent);

        HashMap<String, BamFile> sampleBamMap = new HashMap<>();


        // Loop through each JSON object
        for (int i = 0; i < jsonArray.length(); i++) {

            JSONObject bamEntry = jsonArray.getJSONObject(i);

            seqType es = seqType.valueOf(bamEntry.getString("experimental_strategy"));
            long fs = bamEntry.getLong("file_size");
            String fn = bamEntry.getString("file_name");
            String fid = bamEntry.getString("file_id");
            String s = bamEntry.getString("state");
            String md5 = bamEntry.getString("md5sum");

            String caseId = bamEntry.getJSONArray("associated_entities").getJSONObject(0).getString("case_id");


            BamFile bamFile = new BamFile(es, fs, fn, fid, s, md5);

            //Check if there already is a bam for this sample if so check if this new bam is an improvement.
            if (!sampleBamMap.containsKey(caseId) || sampleBamMap.get(caseId).isOtherBetter(bamFile)) {
                sampleBamMap.put(caseId, bamFile);
            }


        }

        System.out.println("Number of samples with bam: " + sampleBamMap.size());

        int wgsCount = 0;
        int wxsCount = 0;

        for (BamFile bamFile : sampleBamMap.values()) {
            if (bamFile.getEs() == seqType.WGS) {
                wgsCount++;
            } else {
                wxsCount++;
            }
        }

        System.out.println("WGS count: " + wgsCount);
        System.out.println("WXS count: " + wxsCount);


        CSVWriter writer = new CSVWriter(new FileWriter("C:\\Users\\patri\\OneDrive - UMCG\\Documents - Functional Genomics\\Projects\\GutEqtls\\GDC\\metadata.repository.2026-01-26-dna-bam-selected.txt"), '\t', '\0', '\0', "\n");


        String[] outputLine = new String[6];


        for (Map.Entry<String, BamFile> sampleBamEntry : sampleBamMap.entrySet()) {

            BamFile bamFile = sampleBamEntry.getValue();

            int c = 0;
            outputLine[c++] = sampleBamEntry.getKey();
            outputLine[c++] = bamFile.getFid();
            outputLine[c++] = bamFile.getFn();
            outputLine[c++] = bamFile.getMd5();
            outputLine[c++] = String.valueOf(bamFile.getFs());
            outputLine[c++] = bamFile.getS();

            writer.writeNext(outputLine);
        }

        writer.close();

    }


    private enum seqType {
        WXS, WGS;
    }

    private static class BamFile {

        private final seqType es;
        private final long fs;
        private final String fn;
        private final String fid;
        private final String s;
        private final String md5;

        public BamFile(seqType es, long fs, String fn, String fid, String s, String md5) {
            this.es = es;
            this.fs = fs;
            this.fn = fn;
            this.fid = fid;
            this.s = s;
            this.md5 = md5;
        }

        public String getMd5() {
            return md5;
        }

        public seqType getEs() {
            return es;
        }

        public long getFs() {
            return fs;
        }

        public String getFn() {
            return fn;
        }

        public String getFid() {
            return fid;
        }

        public String getS() {
            return s;
        }

        public boolean isOtherBetter(BamFile other) {

            if (this.es == seqType.WXS && other.es == seqType.WGS) {
                return true;
            } else if (this.es == seqType.WGS && other.es == seqType.WXS) {
                return false;
            } else {
                //both have same type so look at file size as proxy for read depth

                return this.fs < other.fs;

            }

        }
    }

}

