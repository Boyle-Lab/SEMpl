/*
 * ContextIntersect.cpp
 *
 *  Created on: Jan 6, 2014
 *      Author: nek3d
 */

#include "ContextIntersect.h"

ContextIntersect::ContextIntersect()
{

}

ContextIntersect::~ContextIntersect()
{

}


bool ContextIntersect::parseCmdArgs(Dataset &data) {

    //Needs to be revised


		if ((data.bedtoolsSettings.options.find("-a")) == string.npos) {
			if (!handle_a(data)) return false;
		}
        if ((data.bedtoolsSettings.options.find("-abam")) == string.npos) {
			if (!handle_abam(data)) return false;
        }
        if ((data.bedtoolsSettings.options.find("-b")) == string.npos) {
			if (!handle_b(data)) return false;
		}
        if ((data.bedtoolsSettings.options.find("-names")) == string.npos) {
			if (!handle_names(data)) return false;
		}
        if ((data.bedtoolsSettings.options.find("-filenames")) == string.npos) {
			if (!handle_filenames(data)) return false;
		}
        if ((data.bedtoolsSettings.options.find("-u")) == string.npos) {
			if (!handle_u()) return false;
        }
        if ((data.bedtoolsSettings.options.find("-f")) == string.npos) {
			if (!handle_f(data)) return false;
        }
        if ((data.bedtoolsSettings.options.find("-F")) == string.npos) {
			if (!handle_F(data)) return false;
        }
        if ((data.bedtoolsSettings.options.find("-wa")) == string.npos) {
			if (!handle_wa()) return false;
        }
        if (( data.bedtoolsSettings.options.find("-wao")) == string.npos) {
			if (!handle_wao(data)) return false;
        }
        if (( data.bedtoolsSettings.options.find("-wb")) == string.npos) {
			if (!handle_wb(data)) return false;
        }
        if (( data.bedtoolsSettings.options.find("-wo")) == string.npos) {
			if (!handle_wo(data)) return false;
        }
        if ((data.bedtoolsSettings.options.find("-c")) == string.npos) {
			if (!handle_c(data)) return false;
        }
        if((data.bedtoolsSettings.options.find("-r")) == string.npos) {
			if (!handle_r(data)) return false;
        }
        if((data.bedtoolsSettings.options.find("-e")) == string.npos) {
			if (!handle_e(data)) return false;
        }
        if ((data.bedtoolsSettings.options.find("-v")) == string.npos) {
			if (!handle_v(data)) return false;
        }
        if ((data.bedtoolsSettings.options.find("-s")) == string.npos) {
			if (!handle_s(data)) return false;
        }
        if ((data.bedtoolsSettings.options.find("-S")) == string.npos) {
			if (!handle_S(data)) return false;
        }
        if ((data.bedtoolsSettings.options.find("-loj")) == string.npos) {
			if (!handle_loj(data)) return false;
        }
	return ContextBase::parseCmdArgs(data);
}



bool ContextIntersect::isValidState()
{
	if (!ContextBase::isValidState()) {
		return false;
	}

	if (_queryFileIdx == -1 || _dbFileIdxs.size() == -0) {
		_errorMsg = "\n***** ERROR: query and database files not specified. *****";
		return false;
	}

	if (getAnyHit() && getNoHit()) {
		_errorMsg = "\n***** ERROR: request either -u for anyHit OR -v for noHit, not both. *****";
		return false;
	}
	if (getWriteCount()) {
		if (getAnyHit()) {
			_errorMsg = "\n***** ERROR: request either -c for writeCount OR -u for anyHit, not both. *****";
			return false;
		}  else if (getWriteB()) {
			_errorMsg = "\n***** ERROR: request either -c for writeCount OR -wb for writeB, not both. *****";
			return false;
		} else if (getQueryFileType() == FileRecordTypeChecker::BAM_FILE_TYPE && !getExplicitBedOutput()) {
			_errorMsg = "\n***** ERROR: writeCount option is not valid with BAM query input, unless bed output is specified with -bed option. *****";
			return false;
		}
	}

	if (_haveFractionA && (_overlapFractionA <= 0.0 || _overlapFractionA > 1.0)) {
		_errorMsg = "\n***** ERROR: -f must be in the range (0.0, 1.0]. *****";
		return false;
	}
	if (_haveFractionB && (_overlapFractionB <= 0.0 || _overlapFractionB > 1.0)) {
		_errorMsg = "\n***** ERROR: -F must be in the range (0.0, 1.0]. *****";
		return false;
	}
	if (_haveFractionB && !_haveFractionA && _eitherFraction) {
		_errorMsg = "\n***** ERROR: -e must be used with -f. *****";
		return false;
	}
	if (_haveFractionB && _reciprocalFraction) {
		_errorMsg = "\n***** ERROR: -r must be used solely with -f. *****";
		return false;
	}
	if (_haveFractionA && _haveFractionB && _reciprocalFraction) {
		_errorMsg = "\n***** ERROR: -r must be used solely with -f. *****";
		return false;
	}

	// if -f and -r are given, then we need to set _overlapFractionB to be the same as _overlapFractionA
	if (_haveFractionA && _reciprocalFraction) {
		setOverlapFractionB(_overlapFractionA);
	}

	if (getUseDBnameTags() && _dbNameTags.size() != _dbFileIdxs.size()) {
		_errorMsg = "\n***** ERROR: Number of database name tags given does not match number of databases. *****";
		return false;

	}


	if (getWriteOverlap()) {

		if (getWriteA()) {
			_errorMsg = "\n***** ERROR: request either -wa for writeA OR -wo for writeOverlap, not both. *****";
			return false;
		} else if (getWriteB()) {
			_errorMsg = "\n***** ERROR: request either -wb for writeB OR -wo for writeOverlap, not both. *****";
			return false;
		}  else if (getWriteCount()) {
			_errorMsg = "\n***** ERROR: request either -c for writeCount OR -wo for writeOverlap, not both. *****";
			return false;
		} else if (getAnyHit()) {
			_errorMsg = "\n***** ERROR: request either -u for anyHit OR -wo for writeOverlap, not both. *****";
			return false;
		} else if (getQueryFileType() == FileRecordTypeChecker::BAM_FILE_TYPE && !getExplicitBedOutput()) {
			_errorMsg = "\n***** ERROR: writeAllOverlap option is not valid with BAM query input, unless bed output is specified with -bed option. *****";
			return false;
		}
	}
	if (getWriteB() || getLeftJoin()) {
		if (getQueryFileType() == FileRecordTypeChecker::BAM_FILE_TYPE && !getExplicitBedOutput()) {
			cerr << endl << "*****" << endl << "*****WARNING: -wb and -loj are ignored with bam input, unless bed output is specified with -bed option." << endl << "*****" << endl;
		}
	}
	if (getSameStrand() && getDiffStrand()) {
		_errorMsg = "\n***** ERROR: request -s for sameStrand, or -S for diffStrand, not both. *****";
		return false;
	}
	if ((getSameStrand() || getDiffStrand()) && !strandedToolSupported()) {
		//error msg set within strandedToolSupported method.
		return false;
	}

	if (getQueryFileType() == FileRecordTypeChecker::BAM_FILE_TYPE && getPrintHeader()) {
		cerr << endl << "*****" << endl << "*****WARNING: -header option is not valid for BAM input." << endl << "*****" << endl;
		setPrintHeader(false);
	}
	if (getAnyHit() || getNoHit() || getWriteCount()) {
		setPrintable(false);
	}

	if (getNoHit() || getWriteCount() || getWriteOverlap()
		|| getWriteAllOverlap() || getLeftJoin())
	{
		setRunToQueryEnd(true);
	}

	if (_files.size()  < 2 ) {
		return false;
	}
	return true;
}

bool ContextIntersect::determineOutputType() {
	if (_outputTypeDetermined) {
		return true;
	}

	//determine the maximum number of database fields.
	for (int i=0; i < getNumDatabaseFiles(); i++) {
		int numFields = getDatabaseFile(i)->getNumFields();
		if ( numFields > _maxNumDatabaseFields) {
			_maxNumDatabaseFields = numFields;
		}
	}

	//If the query is BAM, and bed output wasn't specified, then the output is BAM.
	if (getQueryFileType() == FileRecordTypeChecker::BAM_FILE_TYPE && !getExplicitBedOutput()) {
		setOutputFileType(FileRecordTypeChecker::BAM_FILE_TYPE);
		_outputTypeDetermined = true;
		return true;

	}
	return ContextBase::determineOutputType();

}

bool ContextIntersect::handle_a(Dataset &data)
{
	if (!data.bowtie_output.empty()) {
		_errorMsg = "\n***** ERROR: -a option given, but no query file specified. *****";
		return false;
	}

	//addInputFile("data.bowtie_output");
	_queryFileIdx = getNumInputFiles() -1;
//	markUsed(_i - _skipFirstArgs);
//	_i++;
//	markUsed(_i - _skipFirstArgs);
	return true;
}

bool ContextIntersect::handle_abam(Dataset &data)
{
	if (_argc <= _i+1) {
		_errorMsg = "\n***** ERROR: -abam option given, but no query BAM file specified. *****";
		return false;
	}
	addInputFile(_argv[_i+1]);
	_queryFileIdx = getNumInputFiles() -1;
	markUsed(_i - _skipFirstArgs);
	_i++;
	markUsed(_i - _skipFirstArgs);
	return true;
}

bool ContextIntersect::handle_b(Dataset &data)
{
	if (!data.Dnase.chromosone.empty() ) {
		_errorMsg = "\n***** ERROR: -b option given, but Dnase data not given. *****";
		return false;
	}

//	do {
//		addInputFile(_argv[_i+1]);
//		int fileId = getNumInputFiles() -1;
//		_dbFileIdxs.push_back(fileId);
//		_fileIdsToDbIdxs[fileId] = _dbFileIdxs.size() -1;
//		markUsed(_i - _skipFirstArgs);
//		_i++;
//		markUsed(_i - _skipFirstArgs);
//	} while (_argc > _i+1 && _argv[_i+1][0] != '-');
	return true;
}


bool ContextIntersect::handle_names(Dataset &data)
{
	if (_argc <= _i+1) {
		_errorMsg = "\n***** ERROR: -b option given, but no database names specified. *****";
		return false;
	}

	do {
		addDatabaseNameTag(_argv[_i+1]);
		markUsed(_i - _skipFirstArgs);
		_i++;
		markUsed(_i - _skipFirstArgs);
	} while (_argc > _i+1 && _argv[_i+1][0] != '-');
	setUseDBnameTags(true);
	return true;
}

bool ContextIntersect::handle_filenames(Dataset &data)
{
	markUsed(_i - _skipFirstArgs);
	setUseDBfileNames(true);
	return true;
}

bool ContextIntersect::handle_c(Dataset &data)
{
    setWriteCount(true);
    markUsed(_i - _skipFirstArgs);
    return true;
}


bool ContextIntersect::handle_f(Dataset &data)
{
    if ((_i+1) < _argc) {
        setHaveFractionA(true);
        setOverlapFractionA(atof(_argv[_i + 1]));
        markUsed(_i - _skipFirstArgs);
        _i++;
        markUsed(_i - _skipFirstArgs);
    }
    return true;
}

bool ContextIntersect::handle_F(Dataset &data)
{
    if ((_i+1) < _argc) {
        setHaveFractionB(true);
        setOverlapFractionB(atof(_argv[_i + 1]));
        markUsed(_i - _skipFirstArgs);
        _i++;
        markUsed(_i - _skipFirstArgs);
    }
    return true;
}

bool ContextIntersect::handle_loj(Dataset &data)
{
    setLeftJoin(true);
    markUsed(_i - _skipFirstArgs);
    return true;
}

bool ContextIntersect::handle_r(Dataset &data)
{
	setReciprocalFraction(true);
	markUsed(_i - _skipFirstArgs);
    return true;
}

bool ContextIntersect::handle_e(Dataset &data)
{
	setEitherFraction(true);
	markUsed(_i - _skipFirstArgs);
    return true;
}

bool ContextIntersect::handle_s(Dataset &data)
{
	setSameStrand(true);
	markUsed(_i - _skipFirstArgs);
    return true;
}

bool ContextIntersect::handle_S(Dataset &data)
{
	setDiffStrand(true);
	markUsed(_i - _skipFirstArgs);
    return true;
}

bool ContextIntersect::handle_u()
{
    setAnyHit(true);
    //markUsed(_i - _skipFirstArgs);
    return true;
}

bool ContextIntersect::handle_v(Dataset &data)
{
	setNoHit(true);
	markUsed(_i - _skipFirstArgs);
    return true;
}

bool ContextIntersect::handle_wa()
{
	setWriteA(true);
	//markUsed(_i - _skipFirstArgs);
    return true;
}

bool ContextIntersect::handle_wao(Dataset &data)
{
	setWriteAllOverlap(true);
	setWriteOverlap(true);
	markUsed(_i - _skipFirstArgs);
    return true;
}

bool ContextIntersect::handle_wb(Dataset &data)
{
    setWriteB(true);
     markUsed(_i - _skipFirstArgs);
     return true;
}

bool ContextIntersect::handle_wo(Dataset &data)
{
	setWriteOverlap(true);
	markUsed(_i - _skipFirstArgs);
    return true;
}
