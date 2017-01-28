/*
 File		:	$File: //depot/stillwater-sc/godfa/indexspace/convexhull.go $

 Authors	:	E. Theodore L. Omtzigt
 Date		:	6 April 2014

 Source Control Information:
 Version	:	$Revision: #1 $
 Latest		:	$Date: 2014/04/06 $
 Location	:	$Id: //depot/stillwater-sc/godfa/indexspace/convexhull.go#1 $

 Organization:
		Stillwater Supercomputing, Inc.
		P.O Box 720
		South Freeport, ME 04078-0720

Copyright (c) 2006-2017 E. Theodore L. Omtzigt.  All rights reserved.

Licence      : Stillwater license as defined in this directory

 */
package indexspace

/*
To visualize domains of computation, we need to transform the constraints into the convex hull.
The basic algorithm enumerates all the possible combinations of full rank solutions of subsets of the constraints.
*/