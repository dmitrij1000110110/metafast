package algo;

import it.unimi.dsi.fastutil.longs.LongArrayFIFOQueue;
import org.apache.commons.lang.mutable.MutableLong;
import org.apache.log4j.Logger;
import ru.ifmo.genetics.statistics.Timer;
import ru.ifmo.genetics.structures.map.BigLong2ShortHashMap;
import ru.ifmo.genetics.structures.map.Long2ShortHashMapInterface;
import ru.ifmo.genetics.structures.map.MutableLongShortEntry;
import ru.ifmo.genetics.utils.Misc;
import ru.ifmo.genetics.utils.NumUtils;
import ru.ifmo.genetics.utils.tool.Tool;
import structures.ConnectedComponent;
import structures.ConnectedComponentWithPivots;
import structures.map.BigLong2BitLongaHashMap;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.*;

import static java.lang.Math.max;
import static java.lang.Math.min;

/**
 * Created by -- on 05.02.2021.
 */
public class ComponentsBuilderAroundPivot2 {

    public static List<ConnectedComponentWithPivots> splitStrategy(BigLong2ShortHashMap hm,
                                                                   int k, BigLong2ShortHashMap pivot, BigLong2ShortHashMap pivot2,
                                                                   String statFP, Logger logger, int MIN_PIVOTS, double pivotsQuotient) throws FileNotFoundException {

        logger.debug("split-strategy");

        ComponentsBuilderAroundPivot2 builder = new ComponentsBuilderAroundPivot2(k, statFP, logger, MIN_PIVOTS, pivotsQuotient);
        builder.run(hm, pivot, pivot2);
        return builder.ans;
    }

    final private List<ConnectedComponentWithPivots> ans;
    final int k;
    final String statFP;
    final private Logger logger;

    final int MIN_PIVOTS;

    final double pivotsQuotient;


    public ComponentsBuilderAroundPivot2(int k, String statFP, Logger logger, int MIN_PIVOTS, double pivotsQuotient) {
        this.ans = new ArrayList<ConnectedComponentWithPivots>();
        this.k = k;
        this.statFP = statFP;
        this.logger = logger;
        this.MIN_PIVOTS = MIN_PIVOTS;
        this.pivotsQuotient = pivotsQuotient;
    }

    private void run(BigLong2ShortHashMap hm, BigLong2ShortHashMap pivot, BigLong2ShortHashMap pivot2) throws FileNotFoundException {
        Timer t = new Timer();

        // current component is formed of k-mers with frequency >= 1
        List<ConnectedComponentWithPivots> newComps = findAllComponents(hm, k, pivot, pivot2);

        int ok = 0;


        for (ConnectedComponentWithPivots comp : newComps) {
            ok++;
            ans.add(comp);
        }

        Tool.info(logger, "Found " + NumUtils.groupDigits(ok) + " components");
        Tool.info(logger, "Iteration was finished in " + t);

        Tool.debug(logger, "Memory used: without GC = " + Misc.usedMemoryWithoutRunningGCAsString() + ", " +
                "after it = " + Misc.usedMemoryAsString());

        hm = null;  // for cleaning
        newComps = null;
        Tool.debug(logger, "Memory used after cleaning = " + Misc.usedMemoryAsString() + ", final time = " + t);


        // post processing...
        Tool.debug(logger, "ans.size = " + ans.size());


        Collections.sort(ans);

        PrintWriter statPW = new PrintWriter(statFP);
        statPW.println("# component.no\tcomponent.size\tcomponent.weight\tcomponent.pivotCnt\tcomponent.pivotCnt2\tusedFreqThreshold");
        for (int i = 0; i < ans.size(); i++) {
            ConnectedComponentWithPivots comp = ans.get(i);
            statPW.println((i + 1) + "\t" + comp.size + "\t" + comp.weight + "\t" + comp.pivotCnt + "\t" + comp.pivot2Cnt + "\t" + comp.usedFreqThreshold);
        }
        statPW.close();
    }


    /**
     * Assuming running in one thread for current hm!
     */
    private List<ConnectedComponentWithPivots> findAllComponents(Long2ShortHashMapInterface hm,
                                                                 int k, BigLong2ShortHashMap pivot, BigLong2ShortHashMap pivot2) {
        List<ConnectedComponentWithPivots> ans = new ArrayList<ConnectedComponentWithPivots>();
        LongArrayFIFOQueue queue = new LongArrayFIFOQueue((int) min(1 << 16, hm.size() / 2));
        LongArrayFIFOQueue parent = new LongArrayFIFOQueue((int) min(1 << 16, hm.size() / 2));

        Iterator<MutableLongShortEntry> iterator = pivot.entryIterator();
        while (iterator.hasNext()) {
            MutableLongShortEntry startKmer = iterator.next();
            if (startKmer.getValue() > 0) {    // i.e. if not precessed
                ConnectedComponentWithPivots comp = bfs(hm, startKmer.getKey(), queue, parent, k, pivot, pivot2);
                ans.add(comp);
            }
        }

        return ans;
    }


    /**
     * Breadth-first search to make the traversal of the component.
     * All its kmers are saved to ConnectedComponent.kmers.
     */
    private ConnectedComponentWithPivots bfs(Long2ShortHashMapInterface hm, long startKmer,
                                             LongArrayFIFOQueue queue,
                                             LongArrayFIFOQueue parent, int k, BigLong2ShortHashMap pivot, BigLong2ShortHashMap pivot2) {
        ConnectedComponentWithPivots comp = new ConnectedComponentWithPivots();

        queue.clear();
        parent.clear();

        short value = hm.get(startKmer);
        assert value > 0;
        hm.put(startKmer, (short) -value);  // removing
        comp.add(startKmer, value, 1, 0);
        pivot.put(startKmer, (short) -value);


        // extend to right
        {
            int n_neighbours = 0;
            List<Long> rightNeighbours = new ArrayList<Long>();
            for (long neighbour : KmerOperations.rightNeighbours(startKmer, k)) {
                value = hm.get(neighbour);
                if (value > 0) {
                    rightNeighbours.add(neighbour);
                    n_neighbours++;
                }
            }
            if (n_neighbours == 0) {
                // do nothing
            } else {
                // if single path =>  extend
                if (n_neighbours == 1) {
                    long neighbour = rightNeighbours.get(0);
                    value = hm.get(neighbour);
                    queue.enqueue(neighbour);
                    parent.enqueue(startKmer);
                    hm.put(neighbour, (short) -value);
                    comp.add(neighbour, value, max(0, min(1, pivot.get(neighbour))), max(0, min(1, pivot2.get(neighbour))));
                    if (pivot.get(neighbour) > 0) {
                        pivot.put(neighbour, (short) -value);
                    }
                }
                // if branching path, dfs into each branch to find another pivot or fail
                else {
                    for (long neighbour : rightNeighbours) {
                        List<Long> kmersOnPath = new ArrayList<Long>();
                        boolean goodPath = dfs(neighbour, startKmer, hm, pivot, pivot2, k, kmersOnPath);
                        if (goodPath) {
                            value = hm.get(neighbour);
                            hm.put(neighbour, (short) -value);
                            comp.add(neighbour, value, max(0, min(1, pivot.get(neighbour))), max(0, min(1, pivot2.get(neighbour))));
                            if (pivot.get(neighbour) > 0) {
                                pivot.put(neighbour, (short) -value);
                            }
                            int pathLength = kmersOnPath.size();
                            for (long foundKmer : kmersOnPath) {
                                value = hm.get(foundKmer);
                                comp.add(foundKmer, (short) -value, max(0, min(1, pivot.get(neighbour))), max(0, min(1, pivot2.get(neighbour))));
                            }

                            if (kmersOnPath.size() >= 2) {
                                queue.enqueue(kmersOnPath.get(pathLength - 1));
                                parent.enqueue(kmersOnPath.get(pathLength - 2));
                            } else {
                                if (kmersOnPath.size() == 1) {
                                    queue.enqueue(kmersOnPath.get(pathLength - 1));
                                    parent.enqueue(neighbour);
                                } else {
                                    queue.enqueue(neighbour);
                                    parent.enqueue(startKmer);
                                }
                            }
                        }
                    }
                }
            }
        }

        // extend to left
        {
            int n_neighbours = 0;
            List<Long> leftNeighbours = new ArrayList<Long>();
            for (long neighbour : KmerOperations.leftNeighbours(startKmer, k)) {
                value = hm.get(neighbour);
                if (value > 0) {
                    leftNeighbours.add(neighbour);
                    n_neighbours++;
                }
            }
            if (n_neighbours == 0) {
                // do nothing
            } else {
                // if single path =>  extend
                if (n_neighbours == 1) {
                    long neighbour = leftNeighbours.get(0);
                    value = hm.get(neighbour);
                    queue.enqueue(neighbour);
                    parent.enqueue(startKmer);
                    hm.put(neighbour, (short) -value);
                    comp.add(neighbour, value, max(0, min(1, pivot.get(neighbour))), max(0, min(1, pivot2.get(neighbour))));
                    if (pivot.get(neighbour) > 0) {
                        pivot.put(neighbour, (short) -value);
                    }
                }
                // if branching path, dfs into each branch to find another pivot or fail
                else {
                    for (long neighbour : leftNeighbours) {
                        List<Long> kmersOnPath = new ArrayList<Long>();
                        boolean goodPath = dfs(neighbour, startKmer, hm, pivot, pivot2, k, kmersOnPath);
                        if (goodPath) {
                            value = hm.get(neighbour);
                            hm.put(neighbour, (short) -value);
                            comp.add(neighbour, value, max(0, min(1, pivot.get(neighbour))), max(0, min(1, pivot2.get(neighbour))));
                            if (pivot.get(neighbour) > 0) {
                                pivot.put(neighbour, (short) -value);
                            }
                            int pathLength = kmersOnPath.size();
                            for (long foundKmer : kmersOnPath) {
                                value = hm.get(foundKmer);
                                comp.add(foundKmer, (short) -value, max(0, min(1, pivot.get(neighbour))), max(0, min(1, pivot2.get(neighbour))));
                            }

                            if (kmersOnPath.size() >= 2) {
                                queue.enqueue(kmersOnPath.get(pathLength - 1));
                                parent.enqueue(kmersOnPath.get(pathLength - 2));
                            } else {
                                if (kmersOnPath.size() == 1) {
                                    queue.enqueue(kmersOnPath.get(pathLength - 1));
                                    parent.enqueue(neighbour);
                                } else {
                                    queue.enqueue(neighbour);
                                    parent.enqueue(startKmer);
                                }
                            }
                        }
                    }
                }
            }
        }


        while (queue.size() > 0) {
            long kmer = queue.dequeue();
            long prev = parent.dequeue();

            int right_neighbours = 0;
            List<Long> rightNeighbours = new ArrayList<Long>();
            for (long neighbour : KmerOperations.rightNeighbours(kmer, k)) {
                value = hm.get(neighbour);
                if (value > 0) {
                    rightNeighbours.add(neighbour);
                    right_neighbours++;
                }
            }
            int left_neighbours = 0;
            List<Long> leftNeighbours = new ArrayList<Long>();
            for (long neighbour : KmerOperations.leftNeighbours(kmer, k)) {
                value = hm.get(neighbour);
                if (value > 0) {
                    leftNeighbours.add(neighbour);
                    left_neighbours++;
                }
            }

            int n_neighbours = 0;
            List<Long> neighbours = null;
            for (long val : KmerOperations.leftNeighbours(kmer, k)) {
                if (val == prev) {
                    n_neighbours = right_neighbours;
                    neighbours = rightNeighbours;
                }
            }
            for (long val : KmerOperations.rightNeighbours(kmer, k)) {
                if (val == prev) {
                    n_neighbours = left_neighbours;
                    neighbours = leftNeighbours;
                }
            }

            if (n_neighbours == 0) {
                continue;
                // do nothing
            }
            // if single path =>  extend
            if (n_neighbours == 1) {
                long neighbour = neighbours.get(0);
                value = hm.get(neighbour);
                queue.enqueue(neighbour);
                parent.enqueue(kmer);
                hm.put(neighbour, (short) -value);
                comp.add(neighbour, value, max(0, min(1, pivot.get(neighbour))), max(0, min(1, pivot2.get(neighbour))));
                if (pivot.get(neighbour) > 0) {
                    pivot.put(neighbour, (short) -value);
                }
            }
            // if branching path, dfs into each branch to find another pivot or fail
            else {
                for (long neighbour : neighbours) {
                    List<Long> kmersOnPath = new ArrayList<Long>();
                    boolean goodPath = dfs(neighbour, kmer, hm, pivot, pivot2, k, kmersOnPath);
                    if (goodPath) {
                        value = hm.get(neighbour);
                        hm.put(neighbour, (short) -value);
                        comp.add(neighbour, value, max(0, min(1, pivot.get(neighbour))), max(0, min(1, pivot2.get(neighbour))));
                        if (pivot.get(neighbour) > 0) {
                            pivot.put(neighbour, (short) -value);
                        }
                        int pathLength = kmersOnPath.size();
                        for (long foundKmer : kmersOnPath) {
                            value = hm.get(foundKmer);
                            comp.add(foundKmer, (short) -value, max(0, min(1, pivot.get(neighbour))), max(0, min(1, pivot2.get(neighbour))));
                        }
                        if (kmersOnPath.size() >= 2) {
                            queue.enqueue(kmersOnPath.get(pathLength - 1));
                            parent.enqueue(kmersOnPath.get(pathLength - 2));
                        } else {
                            if (kmersOnPath.size() == 1) {
                                queue.enqueue(kmersOnPath.get(pathLength - 1));
                                parent.enqueue(neighbour);
                            } else {
                                queue.enqueue(neighbour);
                                parent.enqueue(kmer);
                            }
                        }
                    }
                }
            }
        }


        return comp;
    }


    private boolean goodPath(int pivotCnt, int pivot2Cnt) {
        logger.debug(pivotCnt + " " + pivot2Cnt);
        return (double) pivotCnt * pivotsQuotient > (double) pivot2Cnt && pivotCnt >= MIN_PIVOTS;
    }

    private boolean dfs(long startKmer, long parentKmer, Long2ShortHashMapInterface hm,
                        BigLong2ShortHashMap pivot, BigLong2ShortHashMap pivot2, int k, List<Long> kmersOnPath) {
        boolean foundPivot = false;
        int pivotCnt = 0;
        int pivot2Cnt = 0;
        long kmer = startKmer;
        long prev = parentKmer;

        while (true) {
            int right_neighbours = 0;
            List<Long> rightNeighbours = new ArrayList<Long>();
            for (long neighbour : KmerOperations.rightNeighbours(kmer, k)) {
                short value = hm.get(neighbour);
                if (value > 0) {
                    rightNeighbours.add(neighbour);
                    right_neighbours++;
                }
            }
            int left_neighbours = 0;
            List<Long> leftNeighbours = new ArrayList<Long>();
            for (long neighbour : KmerOperations.leftNeighbours(kmer, k)) {
                short value = hm.get(neighbour);
                if (value > 0) {
                    leftNeighbours.add(neighbour);
                    left_neighbours++;
                }
            }

            int n_neighbours = 0;
            List<Long> neighbours = null;
            for (long val : KmerOperations.leftNeighbours(kmer, k)) {
                if (val == prev) {
                    n_neighbours = right_neighbours;
                    neighbours = rightNeighbours;
                }
            }
            for (long val : KmerOperations.rightNeighbours(kmer, k)) {
                if (val == prev) {
                    n_neighbours = left_neighbours;
                    neighbours = leftNeighbours;
                }
            }


            // if single path =>  extend
            if (n_neighbours == 1) {
                long neighbour = neighbours.get(0);
                kmersOnPath.add(neighbour);

                short value = hm.get(neighbour);
                hm.put(neighbour, (short) -value);
                if (pivot.get(neighbour) > 0) {
                    foundPivot = true;
                    pivotCnt++;
                    pivot.put(neighbour, (short) -value);
                }
                if (pivot2.get(neighbour) > 0) {
                    pivot2Cnt++;
                }
                prev = kmer;
                kmer = neighbour;
            }
            // if branching path or no path, stop and return
            else {
                break;
            }

        }

        return goodPath(pivotCnt, pivot2Cnt);
    }


}
