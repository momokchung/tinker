#
#
#  #############################################################
#  ##                                                         ##
#  ##  compile.make  --  compile each of the Tinker routines  ##
#  ##            (Intel Fortran for Linux Version)            ##
#  ##                                                         ##
#  #############################################################
#
#
#  compile all the modules; "sizes" must be first since it is used
#  to set static array dimensions in many of the other modules
#
#
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp sizes.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp action.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp align.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp analyz.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp angang.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp angbnd.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp angpot.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp angtor.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp argue.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp ascii.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp atmlst.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp atomid.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp atoms.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp bath.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp bitor.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp bndpot.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp bndstr.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp bound.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp boxes.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp cell.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp cflux.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp charge.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp chgpen.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp chgpot.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp chgtrn.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp chrono.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp chunks.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp couple.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp ctrpot.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp deriv.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp dipole.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp disgeo.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp disp.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp dma.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp domega.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp dsppot.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp energi.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp ewald.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp expol.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp faces.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp fft.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp fields.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp files.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp fracs.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp freeze.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp gkstuf.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp group.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp hescut.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp hessn.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp hpmf.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp ielscf.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp improp.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp imptor.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp inform.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp inter.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp iounit.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp kanang.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp kangs.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp kantor.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp katoms.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp kbonds.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp kcflux.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp kchrge.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp kcpen.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp kctrn.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp kdipol.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp kdsp.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp kexpl.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp keys.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp khbond.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp kiprop.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp kitors.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp kmulti.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp kopbnd.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp kopdst.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp korbs.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp kpitor.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp kpolr.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp krepl.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp ksolut.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp kstbnd.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp ksttor.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp ktorsn.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp ktrtor.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp kurybr.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp kvdwpr.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp kvdws.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp light.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp limits.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp linmin.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp math.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp mdstuf.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp merck.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp minima.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp molcul.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp moldyn.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp moment.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp mplpot.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp mpole.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp mrecip.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp mutant.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp neigh.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp nonpol.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp nucleo.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp omega.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp opbend.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp opdist.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp openmp.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp orbits.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp output.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp params.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp paths.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp pbstuf.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp pdb.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp phipsi.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp piorbs.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp pistuf.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp pitors.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp pme.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp polar.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp polgrp.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp polopt.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp polpcg.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp polpot.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp poltcg.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp potent.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp potfit.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp ptable.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp qmstuf.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp refer.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp repel.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp reppot.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp resdue.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp restrn.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp rgddyn.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp rigid.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp ring.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp rotbnd.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp rxnfld.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp rxnpot.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp scales.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp sequen.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp shunt.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp socket.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp solpot.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp solute.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp stodyn.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp strbnd.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp strtor.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp syntrn.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp tarray.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp titles.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp torpot.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp tors.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp tortor.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp tree.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp units.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp uprior.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp urey.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp urypot.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp usage.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp valfit.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp vdw.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp vdwpot.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp vibs.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp virial.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp warp.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp xtals.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp zclose.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp zcoord.f
#
#  now compile separately each of the Fortran source files
#
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp active.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp alchemy.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp alterchg.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp alterpol.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp analysis.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp analyze.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp angles.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp anneal.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp arcedit.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp attach.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp baoab.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp bar.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp basefile.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp beeman.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp bicubic.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp bitors.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp bonds.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp born.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp bounds.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp bussi.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp calendar.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp center.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp chkpole.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp chkring.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp chkxyz.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp cholesky.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp clock.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp cluster.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp column.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp command.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp connect.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp connolly.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp control.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp correlate.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp critical.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp crystal.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp cspline.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp cutoffs.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp damping.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp dcflux.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp deflate.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp delete.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp dexpol.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp diagq.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp diffeq.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp diffuse.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp distgeom.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp document.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp dynamic.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp eangang.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp eangang1.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp eangang2.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp eangang3.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp eangle.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp eangle1.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp eangle2.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp eangle3.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp eangtor.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp eangtor1.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp eangtor2.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp eangtor3.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp ebond.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp ebond1.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp ebond2.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp ebond3.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp ebuck.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp ebuck1.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp ebuck2.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp ebuck3.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp echarge.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp echarge1.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp echarge2.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp echarge3.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp echgdpl.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp echgdpl1.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp echgdpl2.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp echgdpl3.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp echgtrn.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp echgtrn1.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp echgtrn2.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp echgtrn3.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp edipole.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp edipole1.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp edipole2.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp edipole3.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp edisp.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp edisp1.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp edisp2.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp edisp3.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp egauss.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp egauss1.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp egauss2.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp egauss3.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp egeom.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp egeom1.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp egeom2.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp egeom3.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp ehal.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp ehal1.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp ehal2.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp ehal3.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp eimprop.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp eimprop1.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp eimprop2.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp eimprop3.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp eimptor.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp eimptor1.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp eimptor2.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp eimptor3.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp elj.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp elj1.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp elj2.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp elj3.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp embed.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp emetal.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp emetal1.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp emetal2.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp emetal3.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp emm3hb.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp emm3hb1.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp emm3hb2.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp emm3hb3.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp empole.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp empole1.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp empole2.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp empole3.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp energy.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp eopbend.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp eopbend1.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp eopbend2.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp eopbend3.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp eopdist.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp eopdist1.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp eopdist2.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp eopdist3.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp epitors.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp epitors1.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp epitors2.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp epitors3.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp epolar.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp epolar1.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp epolar2.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp epolar3.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp erepel.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp erepel1.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp erepel2.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp erepel3.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp erf.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp erxnfld.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp erxnfld1.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp erxnfld2.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp erxnfld3.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp esolv.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp esolv1.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp esolv2.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp esolv3.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp estrbnd.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp estrbnd1.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp estrbnd2.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp estrbnd3.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp estrtor.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp estrtor1.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp estrtor2.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp estrtor3.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp etors.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp etors1.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp etors2.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp etors3.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp etortor.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp etortor1.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp etortor2.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp etortor3.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp eurey.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp eurey1.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp eurey2.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp eurey3.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp evcorr.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp extra.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp extra1.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp extra2.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp extra3.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp fatal.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp fft3d.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp fftpack.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp field.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp final.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp flatten.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp freefix.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp freeunit.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp gda.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp geometry.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp getarc.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp getcart.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp getint.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp getkey.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp getmol.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp getmol2.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp getnumb.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp getpdb.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp getprm.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp getref.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp getstring.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp gettext.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp getword.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp getxyz.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp ghmcstep.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp gradient.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp gradrgd.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp gradrot.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp groups.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp grpline.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp gyrate.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp hessian.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp hessrgd.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp hessrot.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp hybrid.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp image.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp impose.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp induce.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp inertia.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp initatom.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp initial.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp initprm.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp initres.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp initrot.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp insert.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp intedit.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp intxyz.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp invbeta.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp invert.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp jacobi.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp kangang.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp kangle.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp kangtor.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp katom.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp kbond.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp kcharge.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp kchgflx.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp kchgtrn.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp kdipole.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp kdisp.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp kewald.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp kexpol.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp kextra.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp kgeom.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp kimprop.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp kimptor.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp kinetic.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp kmetal.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp kmpole.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp kopbend.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp kopdist.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp korbit.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp kpitors.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp kpolar.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp krepel.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp ksolv.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp kstrbnd.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp kstrtor.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp ktors.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp ktortor.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp kurey.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp kvdw.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp lattice.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp lbfgs.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp lights.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp lusolve.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp makeint.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp makeref.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp makexyz.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp maxwell.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp mdinit.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp mdrest.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp mdsave.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp mdstat.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp mechanic.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp merge.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp minimize.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp minirot.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp minrigid.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp mol2xyz.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp molecule.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp molxyz.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp moments.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp monte.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp mutate.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp nblist.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp newton.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp newtrot.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp nextarg.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp nexttext.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp nose.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp nspline.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp nucleic.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp number.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp numeral.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp numgrad.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp ocvm.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp openend.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp optimize.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp optinit.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp optirot.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp optrigid.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp optsave.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp orbital.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp orient.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp orthog.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp overlap.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp path.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp pdbxyz.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp picalc.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp pmestuf.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp pmpb.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp polarize.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp poledit.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp polymer.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp potential.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp predict.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp pressure.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp prmedit.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp prmkey.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp promo.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp protein.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp prtarc.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp prtdcd.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp prtdyn.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp prterr.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp prtint.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp prtmol2.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp prtpdb.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp prtprm.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp prtseq.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp prtxyz.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp pss.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp pssrigid.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp pssrot.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp qrsolve.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp quatfit.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp radial.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp random.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp rattle.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp readcart.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp readdcd.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp readdyn.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp readgau.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp readgdma.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp readint.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp readmol.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp readmol2.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp readpdb.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp readprm.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp readseq.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp readxyz.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp replica.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp respa.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp rgdstep.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp rings.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp rmsfit.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp rotlist.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp rotpole.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp saddle.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp scan.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp sdstep.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp search.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp server.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp setprm.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp shakeup.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp sigmoid.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp simplex.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp sktstuf.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp sniffer.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp sort.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp spacefill.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp spectrum.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp square.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp suffix.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp superpose.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp surface.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp surfatom.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp switch.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp tcgstuf.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp temper.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp testgrad.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp testhess.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp testpair.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp testpol.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp testrot.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp testvir.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp timer.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp timerot.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp tncg.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp torphase.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp torque.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp torsfit.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp torsions.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp trimtext.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp unitcell.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp valence.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp verlet.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp version.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp vibbig.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp vibrate.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp vibrot.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp volume.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp xtalfit.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp xtalmin.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp xyzatm.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp xyzedit.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp xyzint.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp xyzmol2.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp xyzpdb.f
ifort -c -O3 -xHost -no-ipo -no-prec-div -openmp zatom.f
