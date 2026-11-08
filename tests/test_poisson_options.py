"""Compile the production option validator without the Dune/PDE dependencies."""
import os
from pathlib import Path
import shlex
import subprocess
import tempfile
import unittest


class Options(unittest.TestCase):
    def test_invalid_modes_and_supported_precision(self):
        source=(Path(__file__).resolve().parents[1]/'poisson_bddc_options.hh').read_text()
        validator=source[source.index('struct PoissonBddcOptions'):source.index('inline void configurePoissonTransfer')]
        fixture='''#include <cmath>
#include <stdexcept>
#include <iostream>
namespace Kaskade::BDDC { enum { FACE=4 }; }
'''+validator+'''
void check(bool value) { if (!value) throw std::runtime_error("option validation regression"); }
bool rejected(PoissonBddcOptions o) {
 try { validatePoissonOptions(o); return false; } catch (std::invalid_argument const&) { return true; }
}
int main() {
 PoissonBddcOptions o;
 check(!rejected(o));
 o.mpiRouting=true; check(rejected(o));
 o.mpi=true; check(!rejected(o));
 o.localCompressionBypass=true;check(rejected(o));
 o.compression=true;check(!rejected(o));
 o.graphLiftingBothOrders=true;check(rejected(o));
 o.graphLifting=true;check(!rejected(o));
 o.constantSplit=true;check(rejected(o));
 o.useCg=true;o.relativeTolerance=1e-7;check(!rejected(o));
 for (int bits: {4,6,8,12,16,24,32}) {o.compressionBits=bits;check(!rejected(o));}
 o.compressionBits=48;check(rejected(o));o.compressionBits=16;
 o.flexibleCg=true;check(!rejected(o));o.fgmres=true;check(rejected(o));
 o.flexibleCg=false;check(!rejected(o));
 o.fgmresRestart=0;check(rejected(o));o.fgmresRestart=30;
 o.restrictionPasses=2;check(rejected(o));o.restrictionPasses=1;
 o.relativeTolerance=NAN;check(rejected(o));o.relativeTolerance=1e-7;
 check(!selectSingleNodePoissonBypass(o,2));check(o.compression);
 check(selectSingleNodePoissonBypass(o,1));
 check(!o.compression && !o.localCompressionBypass && !o.constantSplit);
 check(o.originalResidual && o.fgmres && o.mpiRouting);
 check(!rejected(o));
 std::cout << "Poisson mode, precision and routing option guards passed\\n";
}
'''
        with tempfile.TemporaryDirectory() as d:
            cpp=Path(d)/'options.cpp';exe=Path(d)/'options';cpp.write_text(fixture)
            subprocess.run(shlex.split(os.environ.get('CXX','c++'))+['-std=c++20','-Wall','-Wextra','-Werror',str(cpp),'-o',str(exe)],check=True)
            subprocess.run([str(exe)],check=True)


if __name__=='__main__': unittest.main()
