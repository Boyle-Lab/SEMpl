 /*==========================================================================
                SeqAn - The Library for Sequence Analysis
                          http://www.seqan.de
 ============================================================================
  Copyright (C) 2007

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  Lesser General Public License for more details.

 ============================================================================
  $Id: file.h,v 1.4 2009/03/13 14:51:00 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_FILE_H
#define SEQAN_HEADER_FILE_H

//____________________________________________________________________________
// prerequisites

#include <iostream>
#include <climits>
#include <cstdio>
#include <list>
#include <vector>
#include <map>
#include <cmath>

#include "lib/bowtie-1.0.0/SeqAn-1.1/seqan/sequence.h"


//____________________________________________________________________________

#include "lib/bowtie-1.0.0/SeqAn-1.1/seqan/file/file_forwards.h"

#ifdef SEQAN_SWITCH_USE_FORWARDS
#include "lib/bowtie-1.0.0/SeqAn-1.1/seqan/file/file_generated_forwards.h"
#endif

#include "lib/bowtie-1.0.0/SeqAn-1.1/seqan/file/cstream.h"
#include "lib/bowtie-1.0.0/SeqAn-1.1/seqan/file/stream.h"

#include "lib/bowtie-1.0.0/SeqAn-1.1/seqan/file/chunk_collector.h"
#include "lib/bowtie-1.0.0/SeqAn-1.1/seqan/file/meta.h"

//____________________________________________________________________________
// file formats

#include "lib/bowtie-1.0.0/SeqAn-1.1/seqan/file/file_format.h"

#include "lib/bowtie-1.0.0/SeqAn-1.1/seqan/file/stream_algorithms.h"

//file formats for sequences
#include "lib/bowtie-1.0.0/SeqAn-1.1/seqan/file/file_format_raw.h"
#include "lib/bowtie-1.0.0/SeqAn-1.1/seqan/file/file_format_fasta.h"
//#include "lib/bowtie-1.0.0/SeqAn-1.1/seqan/file/file_format_embl.h"
//#include "lib/bowtie-1.0.0/SeqAn-1.1/seqan/file/file_format_genbank.h"

//file formats for alignments
#include "lib/bowtie-1.0.0/SeqAn-1.1/seqan/file/file_format_fasta_align.h"

//others
#include "lib/bowtie-1.0.0/SeqAn-1.1/seqan/file/file_format_cgviz.h"

//____________________________________________________________________________

//#include "lib/bowtie-1.0.0/SeqAn-1.1/seqan/file/file_format_guess.h"

//____________________________________________________________________________
// files

#include "lib/bowtie-1.0.0/SeqAn-1.1/seqan/file/file_base.h"
#include "lib/bowtie-1.0.0/SeqAn-1.1/seqan/file/file_array.h"

//#include "lib/bowtie-1.0.0/SeqAn-1.1/seqan/system.h>	// async file (default file type of File<")
/*#include "lib/bowtie-1.0.0/SeqAn-1.1/seqan/system/file_sync.h"
#include "lib/bowtie-1.0.0/SeqAn-1.1/seqan/system/system_event.h"
#include "lib/bowtie-1.0.0/SeqAn-1.1/seqan/system/file_async.h"
*/
//____________________________________________________________________________
// external string

#endif //#ifndef SEQAN_HEADER_...
