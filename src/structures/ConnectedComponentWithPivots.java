package structures;

import it.unimi.dsi.fastutil.longs.LongArrayList;
import ru.ifmo.genetics.structures.map.BigLong2ShortHashMap;
import ru.ifmo.genetics.utils.tool.ExecutionFailedException;

import java.io.*;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

public class ConnectedComponentWithPivots implements Comparable<ConnectedComponentWithPivots> {

    /**
     * Stores k-mers if the component isn't a big one (less than b2 vertices).
     */
    public List<Long> kmers;

    /**
     * Current component size (number of k-mers)
     */
    public long size;


    public int no;          // id of this component, can be not initialized
    public long weight;     // weight of this component
    public int pivotCnt;
    public int pivot2Cnt;

    /**
     * Threshold used to construct this component!
     * Can be not initialized
     */
    public int usedFreqThreshold;




    /**
     * Stores k-mers with (frequency >= usedFreqThreshold+1) for the following processing.
     */
    public BigLong2ShortHashMap nextHM = null;



    public ConnectedComponentWithPivots() {
        kmers = new LongArrayList();
        size = 0;
        weight = 0;
        pivotCnt = 0;
        pivot2Cnt = 0;
    }

    public void add(long kmer, short w, int cnt, int cnt2) {
        if (cnt2 > 0) {
            //System.out.println(cnt + " " + cnt2);
        }
        kmers.add(kmer);
        size++;
        weight += w;
        pivotCnt += cnt;
        pivot2Cnt += cnt2;
    }
    public void add(long kmer) {
        kmers.add(kmer);
        size++;
    }



    public static void saveComponents(Collection<ConnectedComponentWithPivots> components, String fp) throws IOException {
        DataOutputStream outputStream = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(fp)));
        outputStream.writeInt(components.size());

        for (ConnectedComponentWithPivots component : components) {
            outputStream.writeInt((int) component.size);
            outputStream.writeLong(component.weight);
            for (long kmer : component.kmers) {
                outputStream.writeLong(kmer);
            }
        }

        outputStream.close();
    }

    public static List<ConnectedComponent> loadComponents(File file) throws ExecutionFailedException {
        try {
            DataInputStream inputStream = new DataInputStream(new BufferedInputStream(new FileInputStream(file)));
            int cnt = inputStream.readInt();
            List<ConnectedComponent> res = new ArrayList<ConnectedComponent>(cnt);

            for (int i = 0; i < cnt; i++) {
                int componentSize = inputStream.readInt();
                ConnectedComponent component = new ConnectedComponent();

                component.weight = inputStream.readLong();
                for (int j = 0; j < componentSize; j++) {
                    component.add(inputStream.readLong());
                }
                res.add(component);
                component.no = i + 1;
            }
            inputStream.close();
            return res;
        } catch (FileNotFoundException e) {
            throw new ExecutionFailedException("Can't load components: file not found", e);
        } catch (EOFException e) {
            throw new ExecutionFailedException("Can't load components: file corrupted or format mismatch! " +
                    "Do you set a wrong file?", e);
        } catch (IOException e) {
            throw new ExecutionFailedException("Can't load components: unknown IOException", e);
        }
    }


    @Override
    public int compareTo(ConnectedComponentWithPivots o) {
        int sign = usedFreqThreshold - o.usedFreqThreshold;
        if (sign != 0) {
            return sign;
        }
        sign = -Long.compare(weight, o.weight);
        if (sign != 0) {
            return sign;
        }
        return -Long.compare(size, o.size);
    }
}
