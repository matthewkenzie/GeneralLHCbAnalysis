#!/usr/bin/env python

####################################################
## 																								##
##   Script used for making a generic analysis    ##
##   Will set up the class structure required     ##
## 																								##
####################################################

from optparse import OptionParser
parser = OptionParser()
parser.add_option("-i","--inbranches",default="dat/branches_in.dat",help="Branches to read in. Default: %default")
parser.add_option("-o","--outbranches",default="dat/branches_out.dat",help="Branches to write out. Default: %default")
parser.add_option("-p","--packageLoc",default="../test", help="Location to put the files")
parser.add_option("--firstTime",default=False,action="store_true",help="Use this the first time you run as it will create all the files you need. Be careful about running this later as it may overwrite your hard work!.")
(options,args) = parser.parse_args()

import os
import sys

variableDict = {}

def readBranchesIn():

	f = open(options.inbranches)
	for line in f.readlines():
		if line.startswith('#'): continue
		if line=='' or line=='\n': continue
		els = line.split()
		if len(els)!=3: continue
		vtype = els[0]
		vname = els[1]
		itype = int(els[2])
		variableDict[vname] = [vtype, itype, -1] # -1 means read in

	f.close()

def readBranchesOut():

	f = open(options.outbranches)
	for line in f.readlines():
		if line.startswith('#'): continue
		if line=='' or line=='\n': continue
		els = line.split()
		if len(els)!=3: continue
		vtype = els[0]
		vname = els[1]
		itype = int(els[2])
		if vname in variableDict.keys():
			if vtype!=variableDict[vname][0]:
				sys.exit('ERROR -- Same branch ('+vname+') requested for in and out but with different types!')
			if itype!=variableDict[vname][1]:
				sys.exit('ERROR -- Same branch ('+vname+') requested for in and out but with different types!')
			variableDict[vname] = [vtype, itype, 0] # means read in AND write out
		else:
			variableDict[vname] = [vtype, itype, 1] # 1 means write out

	f.close()

def makePackage():

	if options.firstTime:
		os.system('mkdir -p %s'%options.packageLoc)
		os.system('mkdir -p %s/interface'%options.packageLoc)
		os.system('mkdir -p %s/src'%options.packageLoc)
		os.system('mkdir -p %s/python'%options.packageLoc)
		os.system('mkdir -p %s/dat'%options.packageLoc)
		os.system('cp makefile %s/'%options.packageLoc)

	os.system('touch %s/python/__init__.py'%options.packageLoc)
	os.system('cp %s %s/dat/'%(options.inbranches,options.packageLoc))
	os.system('cp %s %s/dat/'%(options.outbranches,options.packageLoc))

def writeLooperHeader():
	varKeys = variableDict.keys()
	varKeys.sort()

	f = open('%s/interface/Looper.h'%options.packageLoc,'w')
	print 'Writing header to', f.name

	f.write('/////////////////////////////////////\n')
	f.write('//                                 //\n')
	f.write('// Looper.h                        //\n')
	f.write('// Author: Matthew Kenzie          //\n')
	f.write('// Auto-generated                  //\n')
	f.write('// Essentially a wrapper for TTree //\n')
	f.write('//                                 //\n')
	f.write('/////////////////////////////////////\n\n')
	f.write('#ifndef Looper_h\n')
	f.write('#define Looper_h\n')
	f.write('\n#include <iostream>\n')
	f.write('#include <cassert>\n')
	f.write('#include "TChain.h"\n')
	f.write('#include "TTree.h"\n')
	f.write('#include "TBranch.h"\n')
	f.write('#include "TFile.h"\n')
	f.write('#include "TString.h"\n')
	f.write('\n')

	# struct TreeContainer
	f.write('struct TreeContainer {\n')
	f.write('\tTTree *tree;\n')
	f.write('\tTString name;\n')
	f.write('\tLong64_t nentries;\n')
	f.write('\tint itype;\n')
	f.write('\tint sqrts;\n')
	f.write('};\n')
	f.write('\n')

	# class def, cstor, dstor:
	f.write('class Looper {\n')
	f.write('\n\tpublic:\n\n')
	f.write('\t\tLooper(TTree *_outTree, TString _name="Looper");\n')
	f.write('\t\t~Looper();\n\n')

	# class functions
	f.write('\t\t// functions\n')
	f.write('\t\tvoid addTreeContainer(TTree *tree, TString _name, int _itype, int _sqrts);\n')
	f.write('\t\tvoid loadTree(int i);\n')
	f.write('\t\tvoid setOutputBranches();\n')
	f.write('\t\tinline void setFirstEntry(Long64_t ent) { firstEntry = ent; }\n')
	f.write('\t\tinline void setLastEntry(Long64_t ent) { lastEntry = ent; }\n')
	f.write('\t\tLong64_t GetEntries() { return nentries; }\n')
	f.write('\t\tInt_t Fill() { return outTree->Fill(); }\n')

	# class members
	f.write('\n')
	f.write('\t\t// members\n')
	f.write('\t\tTTree *outTree;\n')
	f.write('\t\tTString name;\n')
	f.write('\t\tint itype;\n')
	f.write('\t\tint sqrts;\n')
	f.write('\t\tLong64_t nentries;\n')
	f.write('\t\tLong64_t nbytes;\n')
	f.write('\t\tLong64_t firstEntry;\n')
	f.write('\t\tLong64_t lastEntry;\n')
	f.write('\t\tstd::vector<TreeContainer> treeContainers;\n')
	f.write('\n')

	f.write('\n')
	f.write('\t\t// branch variables and branches\n')
	# write variables
	for key in varKeys:
		var_type = variableDict[key][0]
		ityp = variableDict[key][1]
		line = '\t\t%-15s %s;'%(var_type,key)
		if variableDict[key][2]==1: line += ' // user defined variable'
		f.write(line+'\n')

	# write read in branch pointers
	for key in varKeys:
		var_type = variableDict[key][0]
		ityp = variableDict[key][1]
		line = '\t\t%-15s %s;'%('TBranch','*b_'+key)
		if variableDict[key][2]==1: line += ' // user defined variable'
		f.write(line+'\n')

	f.write('\n};\n\n#endif\n')
	f.close()

def writeLooperSource():
	varKeys = variableDict.keys()
	varKeys.sort()

	f = open('%s/src/Looper.cc'%options.packageLoc,'w')
	print 'Writing source to', f.name

	f.write('/////////////////////////////////////\n')
	f.write('//                                 //\n')
	f.write('// Looper.cc                       //\n')
	f.write('// Author: Matthew Kenzie          //\n')
	f.write('// Auto-generated                  //\n')
	f.write('// Essentially a wrapper for TTree //\n')
	f.write('//                                 //\n')
	f.write('/////////////////////////////////////\n\n')

	f.write('#include "../interface/Looper.h"\n')
	f.write('\nusing namespace std;\n')

	# write constructor implementation
	f.write('\n')
	f.write('Looper::Looper(TTree *_outTree, TString _name):\n')
	f.write('\toutTree(_outTree),\n')
	f.write('\tname(_name),\n')
	f.write('\tnentries(0),\n')
	f.write('\tnbytes(0),\n')
	f.write('\tfirstEntry(-1),\n')
	f.write('\tlastEntry(-1)\n')
	f.write('{\n')
	f.write('\tsetOutputBranches();\n')
	f.write('}\n\n')
	f.write('Looper::~Looper(){}\n\n')

	# addTreeContainer()
	f.write('void Looper::addTreeContainer(TTree *tree, TString _name, int _itype, int _sqrts) {\n\n')
	f.write('\tTreeContainer treeContainer;\n')
	f.write('\ttreeContainer.tree = tree;\n')
	f.write('\ttreeContainer.name = _name;\n')
	f.write('\ttreeContainer.nentries = tree->GetEntries();\n')
	f.write('\ttreeContainer.itype = _itype;\n')
	f.write('\ttreeContainer.sqrts = _sqrts;\n')
	f.write('\n')
	f.write('\ttreeContainers.push_back(treeContainer);\n')
	f.write('\tnentries += treeContainer.nentries;\n')
	f.write('}\n\n')

	# loadTree()
	f.write('void Looper::loadTree(int i) {\n')
	f.write('\tassert(i>=0 && i<treeContainers.size());\n')
	f.write('\n')
	f.write('TTree *tree = treeContainers[i].tree;\n')
	f.write('itype = treeContainers[i].itype;\n')
	f.write('sqrts = treeContainers[i].sqrts;\n')
	f.write('\n')
	for key in varKeys:
		if variableDict[key][2]>0: continue # write only
		line = 'tree->SetBranchAddress("%s", &%s, &b_%s);\n'%(key,key,key)
		if variableDict[key][1]<0: # this is MC only
			line = "if (itype<0) "+line
		if variableDict[key][1]>0: # this is data only
			line = "if (itype>0) "+line
		f.write("\t"+line)
	f.write('}\n\n')

	# setOutputBranches()
	f.write('void Looper::setOutputBranches()\n')
	f.write('{\n')
	f.write('\toutTree->Branch("itype",&itype);\n')
	f.write('\toutTree->Branch("sqrts",&sqrts);\n')
	for key in varKeys:
		if variableDict[key][2]<0: continue # read only
		line = 'outTree->Branch("%s",&%s);\n'%(key,key)
		f.write('\t'+line)
	f.write('}\n\n')

	f.close()

def writeRunnerHeader():
	f = open('%s/interface/Runner.h'%options.packageLoc,'w')
	print 'Writing header to', f.name

	# preamble
	f.write('/////////////////////////////////////\n')
	f.write('//                                 //\n')
	f.write('// Runner.h                        //\n')
	f.write('// Author: Matthew Kenzie          //\n')
	f.write('// Auto-generated                  //\n')
	f.write('// Will run the analysis chain     //\n')
	f.write('//                                 //\n')
	f.write('/////////////////////////////////////\n\n')

	f.write('#ifndef Runner_h\n')
	f.write('#define Runner_h\n')
	f.write('\n')
	f.write('#include <iostream>\n')
	f.write('#include <vector>\n')
	f.write('\n')
	f.write('#include "TString.h"\n')
	f.write('#include "TTree.h"\n')
	f.write('#include "TMath.h"\n')
	f.write('#include "TStopwatch.h"\n\n')
	f.write('#include "../interface/Looper.h"\n')
	f.write('#include "../interface/BaseAnalyser.h"\n')
	f.write('\n')
	# class def
	f.write('class Runner {\n\n')
	f.write('\tpublic:\n\n')
	f.write('\t\tRunner(TTree *_outTree, TString _name="Runner");\n')
	f.write('\t\t~Runner();\n\n')
	f.write('\t\tvoid addLooperTree(TTree *tree, TString name, int itype, int sqrts);\n')
	f.write('\t\tvoid addAnalyser(BaseAnalyser *analyser);\n')
	f.write('\t\tvoid setEntryRange(Long64_t first, Long64_t last) { firstEntry = first; lastEntry = last; } \n')
	f.write('\t\tvoid setFirstEntry(Long64_t ent) { firstEntry = ent; }\n')
	f.write('\t\tvoid setLastEntry(Long64_t ent) { lastEntry = ent; }\n')
	f.write('\t\tvoid printProgressBar(Long64_t jentry, bool isDone=false);\n')
	f.write('\t\tvoid run();\n')
	f.write('\n\tprivate:\n\n')
	f.write('\t\tTString name;\n')
	f.write('\t\tLooper* looper;\n')
	f.write('\t\tstd::vector<BaseAnalyser*> analysers;\n')
	f.write('\t\tLong64_t nentries;\n')
	f.write('\t\tLong64_t firstEntry;\n')
	f.write('\t\tLong64_t lastEntry;\n')
	f.write('\t\tLong64_t naccepted;\n')
	f.write('\t\tTStopwatch timer;\n')
	f.write('\n')
	f.write('};\n\n')
	f.write('#endif\n')

def writeRunnerSource():
	f = open('%s/src/Runner.cc'%options.packageLoc,'w')
	print 'Writing source to', f.name

	f.write('/////////////////////////////////////\n')
	f.write('//                                 //\n')
	f.write('// Runner.cc                       //\n')
	f.write('// Author: Matthew Kenzie          //\n')
	f.write('// Auto-generated                  //\n')
	f.write('// Will run the analysis chain     //\n')
	f.write('//                                 //\n')
	f.write('/////////////////////////////////////\n\n')

	f.write('#include "../interface/Runner.h"\n\n')
	f.write('using namespace std;\n\n')
	# constructor
	f.write('Runner::Runner(TTree *_outTree, TString _name):\n')
	f.write('\tname(_name),\n')
	f.write('\tnentries(0),\n')
	f.write('\tfirstEntry(-1),\n')
	f.write('\tlastEntry(-1),\n')
	f.write('\tnaccepted(0)\n')
	f.write('{\n')
	f.write('\tlooper = new Looper(_outTree,_name);\n')
	f.write('}\n\n')
	#destructor
	f.write('Runner::~Runner(){}\n\n')
	# functions

	# addLooperTree()
	f.write('void Runner::addLooperTree(TTree *tree, TString name, int itype, int sqrts) {\n')
	f.write('\tlooper->addTreeContainer(tree,name,itype,sqrts);\n')
	f.write('\tLong64_t ents = tree->GetEntries();\n')
	f.write('\tnentries += ents;\n')
	f.write('\tcout << Form("%-30s","Runner::addLooperTree()") << " " << "Added LooperTree (" << name.Data() << ") with " << ents << " entries." << endl;\n')
	f.write('}\n\n')

	# addAnalyser()
	f.write('void Runner::addAnalyser(BaseAnalyser *analyser) {\n')
	f.write('\tanalysers.push_back(analyser);\n')
	f.write('\tcout << Form("%-30s","Runner::addAnalyser()") << " " << "Added Analyser (" << analyser->name.Data() << ")." << endl;\n')
	f.write('}\n\n')

	# progress bar
	f.write('void Runner::printProgressBar(Long64_t jentry, bool isDone) {\n')
	f.write('\tdouble percentage = 100.*double(jentry-firstEntry)/double(lastEntry-firstEntry);\n')
	f.write('\tTString prog = "[";\n')
	f.write('\tfor (int i=0; i<=100; i+=2) {\n')
	f.write('\t\tif (percentage>(double(i)-0.001)) prog += "-";\n')
	f.write('\t\telse prog += " ";\n')
	f.write('\t}\n')
	f.write('\tprog += "]";\n\n')
	f.write('\tdouble time = timer.RealTime();\n')
	f.write('\ttimer.Continue();\n')
	f.write('\tdouble timeperevent = time/double(jentry-firstEntry);\n')
	f.write('\tdouble esttimeleft = timeperevent*double(lastEntry-jentry);\n\n')
	f.write('\tif (isDone) percentage = 100.;\n')
	f.write('\tTString summary = Form("%5.1f%% -- %6d/%-6d -- %6.2f ms/ev -- %10.0f secs left",percentage,int(jentry-firstEntry),int(lastEntry-firstEntry),timeperevent*1000.,esttimeleft);\n')
	f.write('\tcout << Form("%-30s","Runner::run()") << " " << prog << " " << summary << "\\r" << flush;\n')
	f.write('}\n\n')

	# run()
	f.write('void Runner::run(){\n\n')
	f.write('\ttimer.Start();\n\n')
	f.write('\tcout << Form("%-30s","Runner::run()") << " " << "Processing Looper (" << looper->name.Data() << ")." << endl;\n')
	f.write('\n')
	f.write('\tfor (unsigned int t=0; t<looper->treeContainers.size(); t++) {\n\n')

	f.write('\t\tLong64_t jentries = looper->treeContainers[t].nentries;\n')
	f.write('\t\tcout << Form("%-30s","Runner::run()") << " " << "Loading tree (" << looper->treeContainers[t].name << ") with entries " << jentries << endl;\n')
	f.write('\t\tlooper->loadTree(t);\n\n')
	f.write('\t\tint cachedFirstEntry = firstEntry;\n')
	f.write('\t\tint cachedLastEntry = lastEntry;\n\n')
	f.write('\t\tif (firstEntry<0) firstEntry=0;\n')
	f.write('\t\tif (lastEntry<0) lastEntry=jentries;\n\n')
	f.write('\t\tfor (Long64_t jentry=0; jentry<jentries; jentry++){\n\n')
	f.write('\t\t\tif (jentry%int(TMath::Ceil(jentries/1000))==0) printProgressBar(jentry);\n')
	f.write('\t\t\tlooper->treeContainers[t].tree->GetEntry(jentry);\n')
	f.write('\t\t\tfor (unsigned a=0; a<analysers.size(); a++){\n')
	f.write('\t\t\t\tif ( ! analysers[a]->AnalyseEvent(*looper) ) continue; // can skip on if the event fails one analysis in the chain\n')
	f.write('\t\t\t\tnaccepted ++;\n')
	f.write('\t\t\t\tlooper->Fill();\n')
	f.write('\t\t\t}\n')
	f.write('\t\t}\n')
	f.write('\t\tprintProgressBar(jentries,true);\n')
	f.write('\t\tcout << endl;\n')
	f.write('\t\tfirstEntry = cachedFirstEntry;\n')
	f.write('\t\tlastEntry = cachedLastEntry;\n')
	f.write('\t}\n')
	f.write('\t\tcout << Form("%-30s","Runner::run()") << " " << "Processing complete. Accepted " << naccepted << " / " << nentries << " events -- " << Form("%6.2f%%",100.*double(naccepted)/double(nentries)) << " efficient" << endl;\n')
	f.write('}\n\n')

def writeAnalyserHeader():
	f = open('%s/interface/BaseAnalyser.h'%options.packageLoc,'w')
	print 'Writing header to', f.name

	f.write('/////////////////////////////////////\n')
	f.write('//                                 //\n')
	f.write('// BaseAnalyser.h                  //\n')
	f.write('// Author: Matthew Kenzie          //\n')
	f.write('// Auto-generated                  //\n')
	f.write('// Will run the analysis chain     //\n')
	f.write('//                                 //\n')
	f.write('/////////////////////////////////////\n\n')

	f.write('#ifndef BaseAnalyser_h\n')
	f.write('#define BaseAnalyser_h\n')
	f.write('\n')
	f.write('#include "TString.h"\n')
	f.write('#include "../interface/Looper.h"\n')
	f.write('\n')
	f.write('class BaseAnalyser {\n')
	f.write('\n\tpublic:\n\n')
	f.write('\t\tBaseAnalyser(TString _name);\n')
	f.write('\t\tvirtual ~BaseAnalyser() = 0;\n')
	f.write('\n')
	f.write('\t\tvirtual bool AnalyseEvent(Looper &l) = 0; // no implementation here\n')
	f.write('\n')
	f.write('\t\tTString name;\n')
	f.write('\n')
	f.write('};\n')
	f.write('\n#endif\n')
	f.close()

def writeAnalyserSource():
	f = open('%s/src/BaseAnalyser.cc'%options.packageLoc,'w')
	print 'Writing source to', f.name

	f.write('/////////////////////////////////////\n')
	f.write('//                                 //\n')
	f.write('// BaseAnalyser.cc                 //\n')
	f.write('// Author: Matthew Kenzie          //\n')
	f.write('// Auto-generated                  //\n')
	f.write('// Will run the analysis chain     //\n')
	f.write('//                                 //\n')
	f.write('/////////////////////////////////////\n\n')

	f.write('#include "../interface/BaseAnalyser.h"\n\n')
	f.write('using namespace std;\n\n')
	f.write('BaseAnalyser::BaseAnalyser(TString _name):\n')
	f.write('\tname(_name)\n')
	f.write('{}\n\n')
	f.write('BaseAnalyser::~BaseAnalyser(){}\n\n')
	f.close()

def writeExampleAnalyserHeader():
	f = open('%s/interface/ExampleAnalyser.h'%options.packageLoc,'w')
	print 'Writing header to', f.name

	f.write('/////////////////////////////////////\n')
	f.write('//                                 //\n')
	f.write('// ExampleAnalyser.h               //\n')
	f.write('// Author: Matthew Kenzie          //\n')
	f.write('// Auto-generated                  //\n')
	f.write('// Will run the analysis chain     //\n')
	f.write('//                                 //\n')
	f.write('/////////////////////////////////////\n\n')

	f.write('#ifndef ExampleAnalyser_h\n')
	f.write('#define ExampleAnalyser_h\n')
	f.write('\n')
	f.write('#include "TString.h"\n')
	f.write('#include "../interface/BaseAnalyser.h"\n')
	f.write('#include "../interface/Looper.h"\n')
	f.write('\n')
	f.write('class ExampleAnalyser : public BaseAnalyser {\n')
	f.write('\n\tpublic:\n\n')
	f.write('\t\tExampleAnalyser(TString _name);\n')
	f.write('\t\t~ExampleAnalyser();\n')
	f.write('\n')
	f.write('\t\tvirtual bool AnalyseEvent(Looper &l); // no implementation here\n')
	f.write('\n')
	f.write('};\n')
	f.write('\n#endif\n')
	f.close()

def writeExampleAnalyserSource():
	f = open('%s/src/ExampleAnalyser.cc'%options.packageLoc,'w')
	print 'Writing source to', f.name

	f.write('/////////////////////////////////////\n')
	f.write('//                                 //\n')
	f.write('// ExampleAnalyser.cc              //\n')
	f.write('// Author: Matthew Kenzie          //\n')
	f.write('// Auto-generated                  //\n')
	f.write('// Will run the analysis chain     //\n')
	f.write('//                                 //\n')
	f.write('/////////////////////////////////////\n\n')

	f.write('#include "../interface/ExampleAnalyser.h"\n\n')
	f.write('using namespace std;\n\n')
	f.write('ExampleAnalyser::ExampleAnalyser(TString _name):\n')
	f.write('\tBaseAnalyser(_name)\n')
	f.write('{}\n\n')
	f.write('ExampleAnalyser::~ExampleAnalyser(){}\n\n')
	f.write('bool ExampleAnalyser::AnalyseEvent(Looper &l){\n')
	f.write('\n\t// do physics here e.g.:\n')
	f.write('\tif ( ! l.B_s0_L0HadronDecision_Dec==1){\n')
	f.write('\t\treturn false;\n')
	f.write('\t}\n\n')
	f.write('\treturn true;\n')
	f.write('}\n')
	f.close()

def writeDatfile():
	f = open('%s/dat/config.dat'%options.packageLoc,'w')
	f.write('####################################################################################################\n') # length is 100 characters
	f.write('##  %-92s  ##\n'%'autogenerated example')
	f.write('##  %-92s  ##\n'%'pass a list of the analyser you want in the first line:')
	f.write('##  %-92s  ##\n'%'   e.g analysers=Analyser,PreSelection,BDTSelection')
	f.write('##  %-92s  ##\n'%'please list all input files in the following format:')
	f.write('##  %-92s  ##\n'%'itype=%itype fname=%fname tname=%tname')
	f.write('##  %-92s  ##\n'%'where:')
	f.write('##  %-92s  ##\n'%'   itype: negative for MC, postive for data with first digit the sqrts')
	f.write('##  %-92s  ##\n'%'          i.e -71 or -72 etc (would be 2011 MC), 81 or 82 etc (would be 2012 data)')
	f.write('##  %-92s  ##\n'%'   fname: path to the input file')
	f.write('##  %-92s  ##\n'%'   tname: name of the input tree')
	f.write('####################################################################################################\n')  # length is 100 characters
	f.write('\n')
	f.write('analysers=ExampleAnalyser\n\n')
	f.write('itype=%d  fname=%-40s  tname=%-20s'%(-71,'exampleinfile.root','exampletreename'))
	f.close()

def writeConfigfile():
	f = open('%s/python/configProducer.py'%options.packageLoc,'w')
	f.write('import ROOT as r\n\n')

	# class infoObj
	f.write('class infoObj:\n')
	f.write('\tdef __init__(self):\n')
	f.write('\t\tself.name = ""\n')
	f.write('\t\tself.fname = ""\n')
	f.write('\t\tself.tname = ""\n')

	# class configProducer
	f.write('\nclass configProducer:\n')
	# constructor
	f.write('\n\tdef __init__(self,runner,datfile,verbose=False):\n\n')
	f.write('\t\tself.runner = runner\n')
	f.write('\t\tself.datfile = datfile\n')
	f.write('\t\tself.verbose = verbose\n')
	f.write('\t\tself.cfgDict = {}\n')
	f.write('\t\tself.analysers = []\n')
	f.write('\n')
	f.write('\t\tself.readDatfile()\n')
	f.write('\t\tif self.verbose: self.printCfg()\n')
	f.write('\t\tself.parseDatfile()\n')
	f.write('\n')

	# getSqrts()
	f.write('\tdef getSqrts(self,itype):\n')
	f.write('\t\tif itype<0: return int(str(itype)[1])\n')
	f.write('\t\telse: return int(str(itype)[0])\n\n')

	# readDatfile()
	f.write('\tdef readDatfile(self):\n\n')
	f.write('\t\tf = open(self.datfile)\n')
	f.write('\t\tfor line in f.readlines():\n')
	f.write('\t\t\tif line.startswith(\'#\'): continue\n')
	f.write('\t\t\tif line==\'\\n\': continue          \n')
	f.write('\t\t\tif len(line.split())==0 and not line.startswith(\'analysers\'): continue\n')
	f.write('\t\t\titype = -999                     \n')
	f.write('\t\t\tline = line.strip(\'\\n\')          \n')
	f.write('\n')
	f.write('\t\t\t# pick up analysers line first\n')
	f.write('\t\t\tif line.startswith(\'analysers\'):\n')
	f.write('\t\t\t\tanalysersLine = line.split(\'=\')[1]\n')
	f.write('\t\t\t\tfor anal in analysersLine.split(\',\'):\n')
	f.write('\t\t\t\t\tself.analysers.append(getattr(r,anal)(anal))\n')
	f.write('\t\t\t\tcontinue\n')
	f.write('\n')
	f.write('\t\t\tels = line.split()               \n')
	f.write('\t\t\t# build cfg dictionary\n')
	f.write('\t\t\tinfo = infoObj()      \n')
	f.write('\t\t\tfor el in els:        \n')
	f.write('\t\t\t\tif el.startswith(\'itype\'):                 \n')
	f.write('\t\t\t\t\titype = int(el.split(\'=\')[1])            \n')
	f.write('\t\t\t\telse:                                      \n')
	f.write('\t\t\t\t\tvarName = el.split(\'=\')[0]               \n')
	f.write('\t\t\t\t\tvarVal = el.split(\'=\')[1]                \n')
	f.write('\n')
	f.write('\t\t\t\t\tsetattr(info,varName,varVal)\n')
	f.write('\n')
	f.write('\t\t\tif itype in self.cfgDict.keys():\n')
	f.write('\t\t\t\tself.cfgDict[itype].append( info )\n')
	f.write('\t\t\telse:\n')
	f.write('\t\t\t\tself.cfgDict[itype] = [ info ] \n')
	f.write('\n')

	# printCfg()
	f.write('\tdef printCfg(self):\n\n')
	f.write('\t\tprint \'%-30s\'%\'configProducer::printCfg()\', \'%-5s  %-5s  %-30s  %-70s  %-20s\'%(\'itype\',\'sqrts\',\'name\',\'fname\',\'tname\')\n')
	f.write('\t\tfor itype, flist in self.cfgDict.items():\n')
	f.write('\t\t\tfor f in flist:\n')
	f.write('\t\t\t\tprint \'%-30s\'%\'\', \'%-5d  %-5d  %-30s  %-70s  %-20s\'%(itype,self.getSqrts(itype),f.name,f.fname, f.tname)\n')
	f.write('\n')
	f.write('\t\tprint \'%-30s\'%\'configProducer::printCfg()\', \'Analysers:\'\n')
	f.write('\t\tfor analyser in self.analysers:\n')
	f.write('\t\t\tprint \'%-30s\'%\'\', analyser.name\n')
	f.write('\n')

	# parseDatfile()
	f.write('\tdef parseDatfile(self):\n\n')
	f.write('\t\tfor itype, flist in self.cfgDict.items():\n')
	f.write('\n')
	f.write('\t\t\tname = flist[0].name\n')
	f.write('\n')
	f.write('\t\t\ttree = r.TChain(name)\n')
	f.write('\n')
	f.write('\t\t\tfor f in flist:\n')
	f.write('\t\t\t\ttree.AddFile(f.fname+\'/\'+f.tname)\n')
	f.write('\t\t\t\tif name != f.name:\n')
	f.write('\t\t\t\t\tsys.exit(\'ERROR -- If the itype is the same for two lines in the datfile, the name must be the same also\')\n')
	f.write('\n')
	f.write('\t\t\tself.runner.addLooperTree(tree,name,itype,self.getSqrts(itype))\n')
	f.write('\n')
	f.write('\t\tfor analyser in self.analysers:\n')
	f.write('\t\t\tself.runner.addAnalyser(analyser)\n')
	f.write('\n')

	f.close()

def writeRunScript():
	f = open('%s/runAnalysis.py'%options.packageLoc,'w')
	f.write('#!/usr/bin/env python\n\n')
	f.write('from optparse import OptionParser\n')
	f.write('parser = OptionParser()\n')
	f.write('parser.add_option("-d","--datfile",default="dat/config.dat",help="Configuration datfile. Default=%default")\n')
	f.write('parser.add_option("-o","--outfile",default="AnalysisOut.root",help="Name of output root file. Default=%default")\n')
	f.write('parser.add_option("-t","--treename",default="AnalysisTree",help="Name of output tree. Default=%default")\n')
	f.write('parser.add_option("-v","--verbose",default=False,action="store_true")\n')
	f.write('(opts,args) = parser.parse_args()\n\n')
	f.write('import ROOT as r\n')
	f.write('r.gSystem.Load("lib/libAnalysis")\n')
	f.write('r.gROOT.SetBatch()\n')
	f.write('print \'%-30s\'%\'runAnalysis.py\', \'Starting analysis\'\n\n')
	f.write('sw = r.TStopwatch()\n')
	f.write('sw.Start()\n\n')
	f.write('from python.configProducer import *\n\n')
	f.write('cfg_file = opts.datfile\n\n')
	f.write('outt = r.TTree(opts.treename,"Analysis Output Tree")\n')
	f.write('\n')
	f.write('runner = r.Runner(outt,"Runner")\n')
	f.write('cfg = configProducer(runner,cfg_file,opts.verbose)\n\n')
	f.write('runner.run()\n')
	f.write('\n')
	f.write('print \'%-30s\'%\'runAnalysis.py\', \'Analysis complete!!!\'\n')
	f.write('print \'%-30s\'%\'runAnalysis.py\', \'Writing tree\', outt.GetName(), \'to output file:\', opts.outfile\n')
	f.write('outf = r.TFile(opts.outfile,"RECREATE")\n')
	f.write('outf.cd()\n')
	f.write('outt.SetDirectory(outf)\n')
	f.write('outt.Write()\n')
	f.write('outf.Close()\n\n')
	f.write('print \'%-30s\'%\'runAnalysis.py\', \'Success!!!\'\n\n')
	f.write('sw.Stop()\n')
	f.write('print \'%-30s\'%\'runAnalysis.py\', \'Took: %4.2f secs (real) %4.2f secs (CPU)\'%(sw.RealTime(),sw.CpuTime())\n')
	f.close()
	os.system('chmod +x %s'%f.name)

## main ##
makePackage()
readBranchesIn()
readBranchesOut()

# headers
writeLooperHeader()

if options.firstTime:
	writeRunnerHeader()
	writeAnalyserHeader()
	writeExampleAnalyserHeader()

# sources
writeLooperSource()

if options.firstTime:
	writeRunnerSource()
	writeAnalyserSource()
	writeExampleAnalyserSource()

# config
if options.firstTime:
	writeConfigfile()
	writeRunScript()
	writeDatfile()
