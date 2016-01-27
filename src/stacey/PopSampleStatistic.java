/*
        Copyright (C) 2015 Graham Jones, www.indriid.com

        This file is part of STACEY.

        STACEY is free software: you can redistribute it and/or modify
        it under the terms of the GNU General Public License as published by
        the Free Software Foundation, either version 3 of the License, or
        (at your option) any later version.

        STACEY is distributed in the hope that it will be useful,
        but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
        GNU General Public License for more details.

        You should have received a copy of the GNU General Public License
        along with STACEY.  If not, see <http://www.gnu.org/licenses/>.
*/

package stacey;

import beast.core.BEASTObject;
import beast.core.Description;
import beast.core.Input;
import beast.core.Loggable;
import beast.core.parameter.RealParameter;
import stacey.util.InverseGammaMixture;

import java.io.PrintStream;

import static beast.util.Randomizer.nextGamma;
import static beast.util.Randomizer.randomChoicePDF;

/**
 *  Created by Graham Jones on 26/09/2014.
 */

@Description("Statistic that samples the the overall population scale factor multiplied by the per-branch inverse gamma mixture.")

public class PopSampleStatistic extends BEASTObject implements Loggable {
    @SuppressWarnings({"CanBeFinal", "WeakerAccess"})
    public Input<RealParameter> popPriorScaleInput =
            new Input<RealParameter>("popPriorScaleInput",
                    "Overall scale for population size",
                    Input.Validate.REQUIRED);


    @SuppressWarnings({"CanBeFinal", "WeakerAccess"})
    public Input<PIOMSCoalescentDistribution> piomscdInput =
            new Input<PIOMSCoalescentDistribution>("piomsCoalDist",
                    "The PIOMSCoalescentDistribution",
                    Input.Validate.REQUIRED);


    public void initAndValidate() {

    }

    @Override
    public void init(PrintStream out) throws Exception {
        out.print("PopSize\t");
    }

    @Override
    public void log(int nSample, PrintStream out) {
        out.print(""+getStatisticValue() + "\t");
    }

    @Override
    public void close(PrintStream out) {
        // nothing to do
    }




    private double getStatisticValue() {
        InverseGammaMixture igm = piomscdInput.get().getInverseGammaMixture();
        double [] a = igm.getAlphas();
        double [] b = igm.getBetas();
        int i = randomChoicePDF(igm.getWeights());
        double x = nextGamma(a[i], b[i]);
        return popPriorScaleInput.get().getValue() / x;
    }

}
