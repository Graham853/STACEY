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

/*
Provides a `union' which represents a subset of the set of species or minimal clusters.
It is a minimal bitset implementation. */


@Description("Low-level utility class for STACEY. <br/>")

public class BitUnion {
private final long bits[];
private final int size;

    // Constructor: makes an empty union, which is a subset of a set with size elements.
    // In application, size is the number of species or minimal clusters.
    public BitUnion(int size) {
        int n = (size+63) / 64;
        bits = new long[n];
        for (int b = 0; b < bits.length; b++) {
            bits[b] = 0;
        }
        this.size = size;
    }

    // for store()
    public void replaceWith(BitUnion x) {
        assert size == x.size;
        System.arraycopy(x.bits, 0, bits, 0, bits.length);
    }

    public void reset() {
        for (int b = 0; b < bits.length; b++) {
            bits[b] = 0;
        }
    }


    // Inserts element indexed by b into a union. No effect if b is already in the union.
    public void insert(int b) {
        assert 0 <= b  &&  b < size;
        bits[b / 64] |= (1L << (b % 64));
    }

    // Returns true if this is contained in x, false if something in this is not in x
    public boolean isContainedIn(BitUnion x) {
        assert size == x.size;
        for (int b = 0; b < bits.length; b++) {
            if ((bits[b] & ~x.bits[b]) != 0) {
                return false;
            }
        }
        return true;
    }


    // Returns true if the intersection of this and x is not empty, false if this and x are disjoint.
    public boolean overlaps(BitUnion x) {
        assert size == x.size;
        for (int b = 0; b < bits.length; b++) {
            if ((bits[b] & x.bits[b]) != 0) {
                return true;
            }
        }
        return false;
    }


    // Replaces this with (this U x)
    public void union(BitUnion x) {
        assert size == x.size;
        for (int b = 0; b < bits.length; b++) {
            bits[b] |= x.bits[b];
        }
    }


    // Returns a textual representation
    public String asText() {
        StringBuilder rep = new StringBuilder();
        rep.append("{");
        for (int b = 0; b < size; b++) {
            String comma = (b==0) ? "" : ",";
            BitUnion x = new BitUnion(size);
            x.insert(b);
            if (x.isContainedIn(this)) {
                rep.append(comma).append(b);
            } else {
                rep.append(comma).append(" ");
                if (b > 9) { rep.append(" "); }
                if (b > 99) { rep.append(" "); }
            }
        }
        rep.append("}");
        return rep.toString();
    }

    public String toString() {
        return asText();
    }


    public int debugNumberBitsSets() {
        int n = 0;
        for (int b = 0; b < size; b++) {
            if ((bits[b / 64] & (1L << (b % 64))) != 0) {
                n++;
            }
        }
        return n;
    }


    public boolean debugBitIsSet(int b) {
        return ((bits[b / 64] & (1L << (b % 64))) != 0);
    }

    public int debugSize() {
        return size;
    }

}
