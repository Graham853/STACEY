package stacey.util;

import beast.evolution.tree.Tree;
import beast.evolution.tree.TreeInterface;
import stacey.GtreeAndCoalFactor;

import java.util.List;

/**
 * Created by Work on 21/12/2015.
 */





// this invariant class just organises the large number of arguments
// passed from PIOMSCoalescentDistribution to ThreadedCalculator
public class PIOMSCInfoForLogPCalc {
    private final FitsHeights fitsHeights;
    private final Bindings bindings;
    private final TreeInterface sTree;
    private final List<GtreeAndCoalFactor> gTreeCFs;
    private final double maxgtreehgt;
    private final boolean [] gTreeFitIsDirty;
    private final boolean [] gTreeFits;
    private final boolean [] gTreeCountIntensityIsDirty;
    private final int [][] coalCounts;
    private final double [][] coalIntensities;
    private final int [] q_j;
    private final double [] gamma_j;
    private final double [] minusLog_r_j;

    public PIOMSCInfoForLogPCalc(FitsHeights fitsHeights, Bindings bindings, TreeInterface sTree, List<GtreeAndCoalFactor> gTreeCFs, double maxgtreehgt,
                                 boolean [] gTreeFitIsDirty, boolean [] gTreeFits, boolean [] gTreeCountIntensityIsDirty,
                                 int [][] coalCounts, double [][] coalIntensities,
                                 int [] q_j, double [] gamma_j, double [] minusLog_r_j) {
        this.fitsHeights = fitsHeights;
        this.bindings = bindings;
        this.sTree = sTree;
        this.gTreeCFs = gTreeCFs;
        this.maxgtreehgt = maxgtreehgt;
        this.gTreeFitIsDirty = gTreeFitIsDirty;
        this.gTreeFits = gTreeFits;
        this.gTreeCountIntensityIsDirty = gTreeCountIntensityIsDirty;
        this.coalCounts = coalCounts;
        this.coalIntensities = coalIntensities;
        this.q_j = q_j;
        this.gamma_j = gamma_j;
        this.minusLog_r_j = minusLog_r_j;
    }

    FitsHeights getFitsHeights() { return fitsHeights; }
    Bindings getBindings() { return bindings; }
    TreeInterface getSTree() { return sTree; }
    List<GtreeAndCoalFactor> getGTreeCFs() { return gTreeCFs; }
    double getMaxgtreehgt() { return maxgtreehgt; }
    boolean [] getGTreeFitIsDirty() { return gTreeFitIsDirty; }
    boolean [] getGTreeFits() { return gTreeFits; }
    boolean [] getGTreeCountIntensityIsDirty() { return gTreeCountIntensityIsDirty; }
    int [][]   getCoalCounts() { return coalCounts; }
    double [][] getCoalIntensities() { return coalIntensities; }
    int [] getQ_j() { return q_j; }
    double [] getGamma_j() { return gamma_j; }
    double [] getMinusLog_r_j() { return minusLog_r_j; }
}
