PROG=calc_mfcc
VER=$(shell grep '^\s*VERSION' ${PROG}.m | sed -e 's/[^0-9.]*\([0-9.]*\)[^0-9.]*/\1/')
DST=${PROG}-v${VER}
#TAR=${DST}.tgz
ZIP=${DST}.zip

ARCH=$(shell ./matlab_arch.sh)
arch=$(shell ./matlab_arch.sh 1)

PROGARCH=${PROG}_${ARCH}
PROGPRJ=${PROG}_prj


MATLAB=/usr/bin/Matlab 
DEPLOYTOOL=/usr/bin/deploytool
#MATLAB=/Applications/MATLAB_R2009a.app/bin/matlab
# MacOS 64 bit (dpwe-macbook)
#DEPLOYTOOL=/Applications/MATLAB_R2010b.app/bin/deploytool
# Linux 64 bit (hog)
#DEPLOYTOOL=/usr/local/MATLAB/R2010b/bin/deploytool 
# Linux 32 bit (cherry)
#DEPLOYTOOL=${MATLAB} -r deploytool

# MCC
#MCC=/Applications/MATLAB_R2010b.app/bin/mcc
MCC=/usr/bin/mcc

# MEX
#MEX=/Applications/MATLAB_R2010b.app/bin/mex
MEX=/usr/bin/mex

DEMOFILE=demo_${PROG}
MAINFILE=${PROG}.m

THUMB=${PROG}_thumb.png

SRCS=${DEMOFILE}.m \
	audspec.m \
	calc_mfcc.m \
        deltas.m \
        fft2melmx.m \
        melfcc.m \
        powspec.m \
        writeasc.m        

DATA=	${THUMB} README \
	test.wav test-ref.out test-htk.out

SUBDIRS= 

EXTRABINS= ReadSound.mexmaci64 ReadSound.mexa64

FORCOMPILE=${PROGPRJ}.prj run_prj_${ARCH}.sh Makefile matlab_arch.sh

DEMOHTML=html/${DEMOFILE}.html
DEMOINDEX=html/index.html

DSTINDEX=${DST}/index.html

# Data to copy to website (referenced from demo_ file) 
# but not otherwise included in packages
NONPKGDATA=\
	${PROG}_MACI64.zip \
	${PROG}_GLNXA64.zip

all: dist

${DEMOHTML}: ${DEMOFILE}.m ${SRCS} ${DATA} 
	${MATLAB} -r "publish ${DEMOFILE}; exit"

${DEMOINDEX}: ${DEMOHTML}
	sed -e 's@<div class="content">@<a href="http://labrosa.ee.columbia.edu/">LabROSA</a> : <a href="http://labrosa.ee.columbia.edu/projects/">Projects</a>: <div class="content"> <IMG SRC="'${THUMB}'" ALIGN="LEFT" HSPACE="10">@' -e 's/amp;auml;/auml;/g' -e 's/@VER@/${VER}/g' < ${DEMOHTML} > ${DEMOINDEX}

compile: ${PROGARCH}.zip

${PROGARCH}.zip: ${SRCS}
	-rm -rf ${PROGPRJ}
	${DEPLOYTOOL} -build ${PROGPRJ}
#	${MCC} -o ${PROGPRJ} -W main:${PROGPRJ} -T link:exe -w enable:specified_file_mismatch -w enable:repeated_file -w enable:switch_ignored -w enable:missing_lib_sentinel -w enable:demo_license -R singleCompThread -R -nodisplay -R -nojvm -v ${MAINFILE}
#	-d ${PROGPRJ}/src 
	mv ${PROGPRJ}/distrib ${PROGPRJ}/${PROGARCH}
	rm ${PROGPRJ}/${PROGARCH}/run_${PROGPRJ}.sh
	cp run_prj_${ARCH}.sh ${PROGPRJ}/${PROGARCH}/run_${PROG}.sh
	rm ${PROGPRJ}/${PROGARCH}/readme.txt
	cp README ${PROGPRJ}/${PROGARCH}/README
	for d in ${SUBDIRS}; do cp -pr $$d ${PROGPRJ}/${PROGARCH}/; done
	for f in ${DATA}; do cp -p $$f ${PROGPRJ}/${PROGARCH}/$$f; done
	cd ${PROGPRJ} && zip -r ${PROGARCH}.zip ${PROGARCH} && cd ..
	mv ${PROGPRJ}/${PROGARCH} ${PROGPRJ}/distrib
	mv ${PROGPRJ}/${PROGARCH}.zip .


test: ${PROGARCH}.zip
	-rm -rf out
	${PROGPRJ}/distrib/run_${PROG}.sh test.wav test.out
	diff test.out test-ref.out

dist: ${SRCS} ${DATA} ${DEMOINDEX} ${EXTRABINS} ${FORCOMPILE}
	rm -rf ${PROG}
	rm -rf ${DST}
	mkdir ${DST}
	cp -pr html/* ${DST}
	rm ${DST}/${DEMOFILE}.html
	for d in ${SUBDIRS}; do cp -pr $$d ${DST}/; done
	for f in ${SRCS} ${DATA} ${EXTRABINS} ${FORCOMPILE}; do cp -p $$f ${DST}/$$f; done
	rm -f ${DST}/*~
	-rm-extended-attribs.sh ${DST}
#	tar cfz ${TAR} ${DST}
	zip -r ${ZIP} ${DST}
# needs to be called PROG (no ver number) not DST on server
	mv ${DST} ${PROG}
	cp -p ${ZIP} ${PROG}
#	cp -p ${PROG}_${ARCH}.zip ${PROG}
	cp -p ${NONPKGDATA} ${PROG}
	scp -pr ${PROG} wool.ee.columbia.edu:public_html/LabROSA/projects/
	scp -pr ${PROG} wool.ee.columbia.edu:wool-public_html/LabROSA/projects/
	scp -pr ${PROG} labrosa.ee.columbia.edu:/var/www/LabROSA/projects/
