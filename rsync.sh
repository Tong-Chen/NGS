#!/bin/bash

while [ 1 == 1 ]
do

	cur=$(date +%Y-%m-%d-%H-%M-%S)
	cd /MPATHD/bak/aliyun500/

	rsync -avzupL --exclude=*.bw --exclude=*.zip --exclude=*.log ct@aliyun500:/var/www/html/devbank4 devbank4 >${cur}.rsync.log 2>&1
	zip -r devbank4.${cur}.zip devbank4 &
	#/bin/rm -rf devbank4.${cur}

	rsync -avzupL --exclude=*.bw --exclude=*.zip --exclude=*.log ct@aliyun500:/var/www/html/M6Adb M6Adb >>${cur}.rsync.log 2>&1
	zip -r M6Adb.${cur}.zip M6Adb &
	#/bin/rm -rf M6Adb.${cur}

	rsync -avzupL --exclude=*.log ct@aliyun500:/var/www/html/IMPObject IMPObject >>${cur}.rsync.log 2>&1
	zip -r IMPObject.${cur}.zip IMPObject &
	#/bin/rm -rf IMPObject.${cur}


	rsync -avzupL ct@aliyun500:/var/www/html/ehbio_doc ehbio_doc >>${cur}.rsync.log 2>&1
	zip -r ehbio_doc.${cur}.zip ehbio_doc
	#/bin/rm -rf ehbio_doc.${cur}

	rsync -avzupL ct@aliyun500:/var/www/html/ImageGP ImageGP >>${cur}.rsync.log 2>&1
	zip -r ImageGP.${cur}.zip ImageGP
	#/bin/rm -rf ImageGP.${cur}

	rsync -avzupL ct@aliyun500:/var/www/html/WheatObject WheatObject >>${cur}.rsync.log 2>&1
	zip -r WheatObject.${cur}.zip WheatObject
	#/bin/rm -rf WheatObject.${cur}

	rsync -avzupL ct@aliyun500:/var/www/html/ehbio_main_page ehbio_main_page >>${cur}.rsync.log 2>&1
	zip -r ehbio_main_page.${cur}.zip ehbio_main_page
	#/bin/rm -rf ehbio_main_page.${cur}

	rsync -avzupL ct@aliyun500:/var/www/html/TCMObject TCMObject >>${cur}.rsync.log 2>&1
	zip -r TCMObject.${cur}.zip TCMObject
	#/bin/rm -rf TCMObject.${cur}

	rsync -avzupL ct@aliyun500:/var/www/html/Training Training >>${cur}.rsync.log 2>&1
	zip -r Training.${cur}.zip Training
	#/bin/rm -rf Training.${cur}

	rsync -avzupL ct@aliyun500:/var/www/html/Train Train >>${cur}.rsync.log 2>&1
	zip -r Train.${cur}.zip Train
	#/bin/rm -rf Train.${cur}

	rsync -avzupL ct@aliyun500:/var/www/html/Esx Esx >>${cur}.rsync.log 2>&1
	zip -r Esx.${cur}.zip Esx
	#/bin/rm -rf Esx.${cur}

	rsync -avzupL ct@aliyun500:/var/www/html/OmicsLab OmicsLab >>${cur}.rsync.log 2>&1
	zip -r OmicsLab.${cur}.zip OmicsLab
	#/bin/rm -rf OmicsLab.${cur}


	sleep 1d

done
