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

package stacey.util;


import beast.core.Description;

import java.util.Arrays;

@Description("Utility class for STACEY which implements a mixture of inverse gamma distributions.")
// It only implements the way the inverse gammas are mixed; it is not a Distribution.
public class InverseGammaMixture {
    private final double [] weights;
    private final double [] alphas;
    private final double [] betas;

    public InverseGammaMixture(int ncs) {
        weights = new double[ncs];
        alphas = new double[ncs];
        betas = new double[ncs];
    }

    public void setWeight(int c, double weight) {
        weights[c] = weight;
    }

    public void setAlpha(int c, double alpha) {
        alphas[c] = alpha;
    }

    public void setBeta(int c, double beta) {
        betas[c] = beta;
    }


    public double [] getWeights() {
        return Arrays.copyOf(weights, weights.length);
    }


    public double[] getAlphas() {
        return Arrays.copyOf(alphas, alphas.length);
    }


    public double[] getBetas() {
        return Arrays.copyOf(betas, betas.length);
    }


    // The mean can be infinite. Finding the mode or median of a mixture is awkward.
    // This returns a value somewhere near the middle of the distribution, but not anything
    // standard.
    public double middlingValue() {
        double mv = 0;
        for (int i = 0; i < weights.length; i++) {
            mv += weights[i] * betas[i] / alphas[i];
        }
        return mv;
    }
}
