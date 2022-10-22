/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package umcg.genetica.io.bin;

import com.mastfrog.util.streams.HashingOutputStream;

import java.io.*;
import java.security.NoSuchAlgorithmException;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

/**
 * @author harm-jan
 */
public class BinaryFile {

    public static final boolean W = true;
    public static final boolean R = false;
    protected final DataOutputStream os;
    protected DataInputStream is;
    protected final File loc;
    protected final boolean writeable;
    private final OutputStream osh;

    public BinaryFile(String loc, boolean mode) throws IOException {
        this(new File(loc), mode, 32 * 1024, false);
    }
    public BinaryFile(File loc, boolean mode) throws IOException {
        this(loc, mode, 32 * 1024, false);
    }

    public BinaryFile(String loc, boolean mode, int buffersize) throws IOException {
        this(new File(loc), mode, buffersize, false);
    }

    public BinaryFile(String loc, boolean mode, int buffersize, boolean useHash) throws IOException {
        this(new File(loc), mode, buffersize, useHash);
    }

    public BinaryFile(File loc, boolean mode, int buffersize, boolean useHash) throws IOException {
        if (loc.getName().trim().length() == 0) {
            throw new IOException("Could not find file: no file specified");
        }
        this.writeable = mode;
        this.loc = loc;

        if (writeable) {
            try {
                is = null;
                if (!useHash) {
                    osh = new FileOutputStream(loc);
                } else {
                    osh = new HashingOutputStream("md5", new FileOutputStream(loc));
                }
                if (loc.getName().endsWith(".gz")) {
                    os = new DataOutputStream(new GZIPOutputStream(new BufferedOutputStream(osh, buffersize)));
                } else {
                    os = new DataOutputStream(new BufferedOutputStream(osh, buffersize));
                }

            } catch (NoSuchAlgorithmException ex) {
                throw new RuntimeException(ex);
            }
        } else {
            if (loc.getName().endsWith(".gz")) {
                is = new DataInputStream(new GZIPInputStream(new BufferedInputStream(new FileInputStream(loc), buffersize)));
            } else {
                is = new DataInputStream(new BufferedInputStream(new FileInputStream(loc), buffersize));
            }

            os = null;
            osh = null;
        }
    }


    public void writeBytes(byte[] v) throws IOException {
        if (writeable) {
            os.write(v);
        } else {
            throw new IOException("File is read only.");
        }
    }

    public void writeInt(int v) throws IOException {
        if (writeable) {
            os.writeInt(v);
        } else {
            throw new IOException("File is read only.");
        }
    }

    public void writeString(String s) throws IOException {
        if (writeable) {
            os.writeUTF(s);
        } else {
            throw new IOException("File is read only.");
        }
    }

    public void writeBool(boolean b) throws IOException {
        if (writeable) {
            os.writeBoolean(b);
        } else {
            throw new IOException("File is read only.");
        }
    }

    public void writeFloat(float f) throws IOException {
        if (writeable) {
            os.writeFloat(f);
        } else {
            throw new IOException("File is read only.");
        }
    }

    public void writeDouble(double d) throws IOException {
        if (writeable) {
            os.writeDouble(d);
        } else {
            throw new IOException("File is read only.");
        }
    }

    public void writeLong(long l) throws IOException {
        if (writeable) {
            os.writeLong(l);
        } else {
            throw new IOException("File is read only.");
        }
    }

    // read functions

    public int readBytes(byte[] b) throws IOException {
        if (writeable) {
            throw new IOException("File is write only.");
        } else {
            return is.read(b);
        }
    }


    public int readInt() throws IOException, EOFException {
        if (writeable) {
            throw new IOException("File is write only.");
        } else {
            return is.readInt();
        }
    }

    public boolean readBool() throws IOException, EOFException {
        if (writeable) {
            throw new IOException("File is write only.");
        } else {
            return is.readBoolean();
        }
    }

    public String readString() throws IOException, EOFException {
        if (writeable) {
            throw new IOException("File is write only.");
        } else {
            return is.readUTF();
        }
    }

    public float readFloat() throws IOException, EOFException {
        if (writeable) {
            throw new IOException("File is write only.");
        } else {
            return is.readFloat();

        }
    }

    public double readDouble() throws IOException, EOFException {
        if (writeable) {
            throw new IOException("File is write only.");
        } else {
            return is.readDouble();
        }
    }

    public long readLong() throws IOException, EOFException {
        if (writeable) {
            throw new IOException("File is write only.");
        } else {
            return is.readLong();
        }
    }

    public void close() throws IOException {
        if (writeable) {
            os.close();
        } else {
            is.close();
        }
    }

    public void writeByte(byte b) throws IOException {
        if (writeable) {
            os.writeByte(b);
        } else {
            throw new IOException("File is read only.");
        }
    }


    public void write(int b) throws IOException {
        os.write(b);
    }

    public byte[] getWrittenHash() throws IOException {
        if (writeable) {
            if (osh instanceof HashingOutputStream) {
                HashingOutputStream hash = (HashingOutputStream) osh;
                return hash.getDigest();
            } else {
                throw new UnsupportedOperationException("OutputStream is not a HashingOutputStream");
            }

        } else {
            return null;
        }
    }

    public int read(byte[] buffer) throws IOException {
        return is.read(buffer);
    }

    public void write(byte[] buffer, int i, int len) throws IOException {
        os.write(buffer, i, len);
    }

    public void write(byte[] alleles0) throws IOException {
        os.write(alleles0);
    }

    public int available() throws IOException {
        return is.available();
    }
}
