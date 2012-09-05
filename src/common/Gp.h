#ifndef GP_H
#define GP_H

#include <iostream>
#include <list>
#include <vector>
#include <stdlib.h>

#ifdef NO_ASSERT
#include <assert.h>
#endif

#ifdef DEBUG
#include <time.h>
#endif

// put in some common header evetually...
#ifndef TRUE
#define TRUE 1
#endif

#ifndef FALSE
#define FALSE 0
#endif

using std::cout;
using std::cerr;
using std::endl;
using std::max;
using std::min;
using std::list;
using std::vector;

#include "MpiStuff.h"
using namespace MpiStuff;

// #############################################

// this is a tommie thig - push it up into Igp...
#define FIRST_PERIODIC_RANGE_BIT 6

class Range
{

private:

  int bits;
  int i_f, i_l;
  int flag;

public:

  double dxyz[3];

  Range(const int i_f, const int i_l, const int bits, const int flag, double * dxyz)
  {
    assert(i_l-i_f+1>=0);
    this->i_f = i_f;
    this->i_l = i_l;
    // bits must be zero or positive?...
    // I guess not (why not use the sign bit as well), but for now...
    assert(bits>=0);
    this->bits = bits;

    // flag is used to set periodicity...
    this->flag = flag;

    // transform...
    if (dxyz==NULL)
    {
      for (int i = 0; i<3; i++)
        this->dxyz[i] = 0.0;
    }
    else
    {
      for (int i = 0; i<3; i++)
        this->dxyz[i] = dxyz[i];
    }
  }

  Range(const int i_f, const int i_l, const int bits)
  {
    assert(i_l-i_f+1>=0);
    this->i_f = i_f;
    this->i_l = i_l;
    // bits must be zero or positive?...
    // I guess not (why not use the sign bit as well), but for now...
    assert(bits>=0);
    this->bits = bits;
    // flag gets zero's by default...
    flag = 0;
  }

  int getBits() const
  {
    return (bits);
  }
  int getFlag() const
  {
    return (flag);
  }
  int getIndexFirst() const
  {
    return (i_f);
  }
  int getIndexLast() const
  {
    return (i_l);
  }
  void setIndexFirst(const int i_f)
  {
    this->i_f = i_f;
  }
  void setIndexLast(const int i_l)
  {
    this->i_l = i_l;
  }

  // for the Tommie heirarchy...
  // the specific periodic or other transform associated with a
  // given range is indicated by the bits. if there is no transform,
  // then bits == 0, otherwise transforms and their inverse transforms
  // use bit pairs 0:1, 2:3, 4:5, etc. To support the MAX_NO_PERIODIC_DATA
  // routine and treat FULL paired boundaries of Igp differently, formally
  // periodic boundaries do not start with bit pair 0:1, but rather 6:7,
  // which is defined as FIRST_PERIODIC_RANGE_BIT above.
  //
  // for the Ugp hierarchy, there is no full, and any non-zero bit represents
  // a periodic transform: ie: isPeriodic() == bits > 0.

  //int isPeriodicTommie() const { return( bits >= (1<<FIRST_PERIODIC_RANGE_BIT) ); }

  // any bit indicates periodicity in UGP...
  int isPeriodic() const
  {
    return (bits>0);
  }
  int isInternal() const
  {
    return (bits==0);
  }

};

class Prcomm
{

private:

  int nbrRank;

public:

  list<Range> packRangeList, unpackRangeList;
  vector<int> packIndexVec, unpackIndexVec;

  int packBufferIntSize, unpackBufferIntSize, packBufferDoubleSize, unpackBufferDoubleSize;
  int * packBufferInt;
  int * unpackBufferInt;
  double * packBufferDouble;
  double * unpackBufferDouble;

  // for variable-sized pack/unpack...
  int npack_v, nunpack_v;

  // for shared memory copying...
  double * packBufferDoublePtr;

  // use these to store changes
  vector<int> packIndexVec2, unpackIndexVec2;

  Prcomm(int rank)
  {

    //cout << mpi_rank << ": new Prcomm for nbrRank: " << rank << endl;

    nbrRank = rank;

    packBufferInt = NULL;
    unpackBufferInt = NULL;
    packBufferDouble = NULL;
    unpackBufferDouble = NULL;

    packBufferIntSize = 0;
    unpackBufferIntSize = 0;
    packBufferDoubleSize = 0;
    unpackBufferDoubleSize = 0;

    packBufferDoublePtr = NULL;

  }

  int getNbrRank() const
  {
    return (nbrRank);
  }

  // there is some convergence neccessary on how exactly ranges are built and their
  // periodic transforms supported. I am trying to make the Range carry the type
  // of the transform along with the transform data, and then the Ugp parent
  // defines a virtual translation and rotation for each flag of the range.

  Range * addPackRange(const int i_f, const int i_l, const int bits, const int flag, double * dxyz)
  {

    // an internal pack range has zero bits - there should only be one,
    // so make sure a zero range does not already exist...

    list<Range>::iterator ri;
    for (ri = packRangeList.begin(); ri!=packRangeList.end(); ri++)
    {
      if (ri->getBits()==bits)
      {
        cerr<<"Error: packRange with bits already present: "<<bits<<endl;
        throw(-1);
      }
    }

    // also, if the bits are zero (indicating no transform), this should
    // be the first pack range...
    if (bits==0)
    {
      assert(packRangeList.size()==0);
    }

    if (dxyz==NULL)
    {
      assert(bits==0);
    }

    // add it to the end, returning a pointer
    packRangeList.push_back(Range(i_f, i_l, bits, flag, dxyz));
    return (&(packRangeList.back()));

  }

  Range * addUnpackRange(const int i_f, const int i_l, const int bits, const int flag, double * dxyz)
  {

    // an internal unpack range has zero bits - there should only be one,
    // so make sure a zero range does not already exist...

    list<Range>::iterator ri;
    for (ri = unpackRangeList.begin(); ri!=unpackRangeList.end(); ri++)
    {
      if (ri->getBits()==bits)
      {
        cerr<<"Error: unpackRange with bits already present: "<<bits<<endl;
        throw(-1);
      }
    }

    // also, if the bits are zero (indicating no transform), this should
    // be the first unpack range...
    if (bits==0)
    {
      assert(unpackRangeList.size()==0);
    }

    if (dxyz==NULL)
    {
      assert(bits==0);
    }

    // add it to the end, returning a pointer
    unpackRangeList.push_back(Range(i_f, i_l, bits, flag, dxyz));
    return (&(unpackRangeList.back()));

  }

  Range * getPackRange(int bits)
  {

    // every range is uniquely determined by its bits. Pack Ranges
    // are also stored such that their flags are in sorted order - i.e.
    // this normally means internal boundaries first, then periodic
    // boundaries. Note that this has implications on how pack
    // and unpack ranges are flagged for periodic boundaries.

    int i_f = packIndexVec.size();

    list<Range>::iterator ri;
    for (ri = packRangeList.begin(); ri!=packRangeList.end(); ri++)
    {
      if (ri->getBits()==bits)
      {
        return (&(*ri));
      }
      else if (ri->getBits()>bits)
      {
        i_f = ri->getIndexFirst();
        break;
      }
    }

    // instantiate a new 0-length pack range starting at i_f, which
    // is set to the start of the next range or this prcomm's packIndexVec...

    packRangeList.insert(ri, Range(i_f, i_f-1, bits)); // inserts BEFORE ri
    ri--; // backs up to the new element
    return (&(*ri));

  }

  Range * getUnpackRange(int bits)
  {

    // see note above. Unpack ranges have a slightly different sorting order in
    // that range 0 is first, then range 2, range 1, range 4 range 3, etc...
    // i.e. odd/even pairs are switched. This supports proper pack/unpack
    // ordering.

    int i_f = unpackIndexVec.size();

    list<Range>::iterator ri;
    for (ri = unpackRangeList.begin(); ri!=unpackRangeList.end(); ri++)
    {
      if (ri->getBits()==bits)
      {
        // a match...
        return (&(*ri));
      }
      else if (bits==0)
      {
        // new one has to be the first one, so if
        // we did not match above, then add the range with
        // zer bits first...
        i_f = ri->getIndexFirst();
        break;
      }
      else if (ri->getBits()!=0)
      {
        // compare which bits are set...
        // for this routine we expect one and ONLY one range to be set...
        int bit = 0;
        int bit_r = 0;
        int ii;
        for (ii = 0; ii<32; ii++)
        {
          if (bits==(1<<ii)) bit = ii;
          if (ri->getBits()==(1<<ii)) bit_r = ii;
        }
        assert(bit>=0);
        assert(bit_r>=0);
        assert(bit!=bit_r); // done above
        // if the new bit is odd, it goes BEFORE the smaller even bit
        // i.e. 1:0,3:2,5:4, etc...
        if (((bit%2!=0)&&(bit<=bit_r+1))||((bit%2==0)&&(bit<bit_r-1)))
        {
          i_f = ri->getIndexFirst();
          break;
        }
      }
    }

    // instantiate a new 0-length unpack range starting at i_f, which
    // is set to the start of the next range or this prcomm's unpackIndexVec...

    unpackRangeList.insert(ri, Range(i_f, i_f-1, bits)); // inserts BEFORE ri
    ri--; // backs up to the new element
    return (&(*ri));

  }

  void addElementsToRangeEnd(int * elements, const int n, Range * range)
  {

    // add elements to the end of the passed range - could be pack or unpack...

    list<Range>::iterator ri;
    for (ri = packRangeList.begin(); ri!=packRangeList.end(); ri++)
    {
      if (&(*ri)==range)
      {
        // make space in the packIndexVec for these "n" elements...
        int n_old = packIndexVec.size();
        packIndexVec.resize(n_old+n);
        // we will insert at the end of this range, so everything
        // above must be shifted up by "n"...
        int i_l = range->getIndexLast();
        int i;
        for (i = n_old-1; i>i_l; i--)
          packIndexVec[i+n] = packIndexVec[i];
        for (i = 0; i<n; i++)
          packIndexVec[i_l+1+i] = elements[i];
        range->setIndexLast(i_l+n);
        // and increment the remaining range indices...
        while (++ri!=packRangeList.end())
        {
          ri->setIndexFirst(ri->getIndexFirst()+n);
          ri->setIndexLast(ri->getIndexLast()+n);
        }
        return;
      }
    }

    for (ri = unpackRangeList.begin(); ri!=unpackRangeList.end(); ri++)
    {
      if (&(*ri)==range)
      {
        // make space in the unpackIndexVec for these "n" elements...
        int n_old = unpackIndexVec.size();
        unpackIndexVec.resize(n_old+n);
        // we will insert at the end of this range, so everything
        // above must be shifted up by "n"...
        int i_l = range->getIndexLast();
        int i;
        for (i = n_old-1; i>i_l; i--)
          unpackIndexVec[i+n] = unpackIndexVec[i];
        for (i = 0; i<n; i++)
          unpackIndexVec[i_l+1+i] = elements[i];
        range->setIndexLast(i_l+n);
        // and increment the remaining range indices...
        while (++ri!=unpackRangeList.end())
        {
          ri->setIndexFirst(ri->getIndexFirst()+n);
          ri->setIndexLast(ri->getIndexLast()+n);
        }
        return;
      }
    }

    cerr<<"Error: range not found in prcomm 1."<<endl;
    throw(-1);

  }

  void addElementsToRangeStart(int * elements, const int n, Range * range)
  {

    // add elements to the start of the passed range...

    list<Range>::iterator ri;
    for (ri = packRangeList.begin(); ri!=packRangeList.end(); ri++)
    {
      if (&(*ri)==range)
      {
        // make space in the packIndexVec for these "n" elements...
        int n_old = packIndexVec.size();
        packIndexVec.resize(n_old+n);
        // we will insert at the end of this range, so everything
        // above must be shifted up by "n"...
        int i_f = range->getIndexFirst();
        int i;
        for (i = n_old-1; i>=i_f; i--)
          packIndexVec[i+n] = packIndexVec[i];
        for (i = 0; i<n; i++)
          packIndexVec[i_f+i] = elements[i];
        range->setIndexLast(range->getIndexLast()+n);
        // and increment the remaining range indices...
        while (++ri!=packRangeList.end())
        {
          ri->setIndexFirst(ri->getIndexFirst()+n);
          ri->setIndexLast(ri->getIndexLast()+n);
        }
        return;
      }
    }

    for (ri = unpackRangeList.begin(); ri!=unpackRangeList.end(); ri++)
    {
      if (&(*ri)==range)
      {
        // make space in the unpackIndexVec for these "n" elements...
        int n_old = unpackIndexVec.size();
        unpackIndexVec.resize(n_old+n);
        // we will insert at the end of this range, so everything
        // above must be shifted up by "n"...
        int i_f = range->getIndexFirst();
        int i;
        for (i = n_old-1; i>=i_f; i--)
          unpackIndexVec[i+n] = unpackIndexVec[i];
        for (i = 0; i<n; i++)
          unpackIndexVec[i_f+i] = elements[i];
        range->setIndexLast(range->getIndexLast()+n);
        // and increment the remaining range indices...
        while (++ri!=unpackRangeList.end())
        {
          ri->setIndexFirst(ri->getIndexFirst()+n);
          ri->setIndexLast(ri->getIndexLast()+n);
        }
        return;
      }
    }

    cerr<<"Error: range not found in prcomm 2."<<endl;
    throw(-1);

  }

  void ensurePackBufferIntSize(const int npack)
  {

    if (npack>packBufferIntSize)
    {
      int * ptr;
      if ((ptr = (int*) realloc(packBufferInt, sizeof(int)*npack))==NULL)
      {
        cerr<<"Error: realloc failed in ensurePackBufferIntSize: "<<npack<<" "<<packBufferIntSize<<endl;
        throw(-1);
      }
      packBufferInt = ptr;
      packBufferIntSize = npack;
    }

  }

  void ensurePackBufferDoubleSize(const int npack)
  {

    if (npack>packBufferDoubleSize)
    {
      double * ptr;
      if ((ptr = (double*) realloc(packBufferDouble, sizeof(double)*npack))==NULL)
      {
        cerr<<"Error: realloc failed in ensurePackBufferDoubleSize: "<<npack<<" "<<packBufferDoubleSize<<endl;
        throw(-1);
      }
      packBufferDouble = ptr;
      packBufferDoubleSize = npack;
    }

  }

  void ensureUnpackBufferIntSize(const int nunpack)
  {

    if (nunpack>unpackBufferIntSize)
    {
      int * ptr;
      if ((ptr = (int*) realloc(unpackBufferInt, sizeof(int)*nunpack))==NULL)
      {
        cerr<<"Error: realloc failed in ensureUnpackBufferIntSize: "<<nunpack<<" "<<unpackBufferIntSize<<endl;
        throw(-1);
      }
      unpackBufferInt = ptr;
      unpackBufferIntSize = nunpack;
    }

  }

  void ensureUnpackBufferDoubleSize(const int nunpack)
  {

    if (nunpack>unpackBufferDoubleSize)
    {
      double * ptr;
      if ((ptr = (double*) realloc(unpackBufferDouble, sizeof(double)*nunpack))==NULL)
      {
        cerr<<"Error: realloc failed in ensureUnpackBufferDoubleSize: "<<nunpack<<" "<<unpackBufferDoubleSize<<endl;
        throw(-1);
      }
      unpackBufferDouble = ptr;
      unpackBufferDoubleSize = nunpack;
    }

  }

};

class PrcommListWrapper
{

  // this class wraps the prcomm lists and adds data neccessary
  // to support non-blocking messaging...

public:

  list<Prcomm> * prcommListPtr;

  int nonSelfPeriodicSize;

  int action;
  void * ptr;

  MPI_Request * sendRequestArray;
  MPI_Request * recvRequestArray;
  MPI_Status * statusArray;
  vector<Prcomm*> prcommPtrVec;

  PrcommListWrapper()
  {

    prcommListPtr = NULL;

    nonSelfPeriodicSize = -1;

    action = -1;
    ptr = NULL;

    sendRequestArray = NULL;
    recvRequestArray = NULL;
    statusArray = NULL;
  }

  ~PrcommListWrapper()
  {
    if (sendRequestArray!=NULL) delete[] sendRequestArray;
    if (recvRequestArray!=NULL) delete[] recvRequestArray;
  }

  void setPrcommList(list<Prcomm>& prcommList)
  {

    assert(prcommListPtr==NULL);
    prcommListPtr = &prcommList;

  }

};

// update actions...

#define REPLACE_DATA             1
#define ADD_DATA                 2
#define MIN_DATA                 3
#define MAX_DATA                 4
#define REPLACE_TRANSLATE_DATA   5
#define REPLACE_ROTATE_DATA      6
#define ADD_ROTATE_DATA          7
#define BITWISE_OR_DATA          8
#define MIN_NO_PERIODIC_DATA     9
#define MAX_NO_PERIODIC_DATA     10
#define ADD_NO_PERIODIC_DATA     11
#define ADD_TRANSLATE_DATA       12
#define SUBTRACT_DATA            13
#define SUBTRACT_ROTATE_DATA     14
#define NO_CHANGE_DATA           15

class Gp
{

protected:

  list<Prcomm> facePrcommList;
  list<Prcomm> edgePrcommList;
  list<Prcomm> nodePrcommList;

  PrcommListWrapper facePrcommListWrapper;
  PrcommListWrapper edgePrcommListWrapper;
  PrcommListWrapper nodePrcommListWrapper;

public:

  Gp()
  {

    facePrcommListWrapper.setPrcommList(facePrcommList);
    edgePrcommListWrapper.setPrcommList(edgePrcommList);
    nodePrcommListWrapper.setPrcommList(nodePrcommList);

  }

  Prcomm * getFacePrcomm(int rank)
  {
    list<Prcomm>::iterator prcomm = facePrcommList.begin();
    while (prcomm!=facePrcommList.end())
    {
      if (prcomm->getNbrRank()==rank) return (&(*prcomm));
      prcomm++;
    }
    facePrcommList.push_back(Prcomm(rank));
    return (&(facePrcommList.back()));
  }

  Prcomm * getEdgePrcomm(int rank)
  {
    list<Prcomm>::iterator prcomm = edgePrcommList.begin();
    while (prcomm!=edgePrcommList.end())
    {
      if (prcomm->getNbrRank()==rank) return (&(*prcomm));
      prcomm++;
    }
    edgePrcommList.push_back(Prcomm(rank));
    return (&(edgePrcommList.back()));
  }

  Prcomm * getNodePrcomm(int rank)
  {
    list<Prcomm>::iterator prcomm = nodePrcommList.begin();
    while (prcomm!=nodePrcommList.end())
    {
      if (prcomm->getNbrRank()==rank) return (&(*prcomm));
      prcomm++;
    }
    nodePrcommList.push_back(Prcomm(rank));
    return (&(nodePrcommList.back()));
  }

protected:

  // when you inherit Gp, you must provide these routines to
  // tell the code how to handle periodicity...

  virtual void rangeTranslateVector(double * xi, const int n, const int flag, double * transform) = 0;
  virtual void rangeRotateVector(double * vi, const int n, const int flag, double * transform) = 0;
  virtual void rangeRotateTensor(double * tij, const int n, const int flag, double * transform) = 0;
  virtual void rangeRotateSymmetricTensor(double * tij, const int n, const int flag, double * transform) = 0;

  // for faces

  void updateFaData(int * d, const int action)
  {
    updateI1(d, action, facePrcommList);
  }

  void updateFaI1(int * d, const int action)
  {
    updateI1(d, action, facePrcommList);
  }

  void updateFaI2(int(*d)[3], const int action)
  {
    updateI2(d, action, facePrcommList);
  }

  void updateFaData(double * d, const int action)
  {
    updateR1(d, action, facePrcommList);
  }

  void updateFaR1(double * d, const int action)
  {
    updateR1(d, action, facePrcommList);
  }

  void updateFaData(double(*d)[3], const int action)
  {
    updateR2(d, action, facePrcommList);
  }

  void updateFaR2(double(*d)[3], const int action)
  {
    updateR2(d, action, facePrcommList);
  }

  void updateFaR3(double(*d)[3][3], const int action)
  {
    updateR3(d, action, facePrcommList);
  }

  void updateFaR2Reverse(double(*d)[3], const int action)
  {
    updateR2Reverse(d, action, facePrcommList);
  }

  void updateFaSymmetricR3(double(*diag)[3], double(*offd)[3], const int action)
  {
    updateSymmetricR3(diag, offd, action, facePrcommList);
  }

  // for nodes

  void updateNoData(int * d, const int action)
  {
    updateI1(d, action, nodePrcommList);
  }

  void updateNoI1(int * d, const int action)
  {
    updateI1(d, action, nodePrcommList);
  }

  void updateNoI1Reverse(int * d, const int action)
  {
    updateI1Reverse(d, action, nodePrcommList);
  }

  void updateNoI2(int(*d)[3], const int action)
  {
    updateI2(d, action, nodePrcommList);
  }

  void updateNoI2Reverse(int(*d)[3], const int action)
  {
    updateI2Reverse(d, action, nodePrcommList);
  }

  void updateNoData(double * d, const int action)
  {
    updateR1(d, action, nodePrcommList);
  }

  void updateNoR1(double * d, const int action)
  {
    updateR1(d, action, nodePrcommList);
  }

  void updateNoData(double(*d)[3], const int action)
  {
    // make sure the nodePrcommList buffers are not currently
    // in use by another non-blocking message exchange...
    assert(nodePrcommListWrapper.ptr==NULL);
    updateR2(d, action, nodePrcommList);
  }

  void updateNoR2(double(*d)[3], const int action)
  {
    // make sure the nodePrcommList buffers are not currently
    // in use by another non-blocking message exchange...
    assert(nodePrcommListWrapper.ptr==NULL);
    updateR2(d, action, nodePrcommList);
  }

  void updateNoR1Start(double *d, const int action)
  {

    // make sure the nodePrcommList buffers are not currently
    // in use by another non-blocking message exchange...
    assert(nodePrcommListWrapper.ptr==NULL);

    // set the pointer and action for consistency checking in Finish
    nodePrcommListWrapper.ptr = (void*) d;
    nodePrcommListWrapper.action = action;

    // do the packing and post the mpi messages...
    updateR1Start(d, action, nodePrcommListWrapper);

  }

  void updateNoR1Finish(double *d, const int action)
  {

    // consistency check...
    assert(nodePrcommListWrapper.ptr==(void*) d);
    assert(nodePrcommListWrapper.action==action);

    // do the work...
    updateR1Finish(d, action, nodePrcommListWrapper);

    // and reset the pointer...
    nodePrcommListWrapper.ptr = NULL;

  }

  void updateNoR2Start(double(*d)[3], const int action)
  {

    // make sure the nodePrcommList buffers are not currently
    // in use by another non-blocking message exchange...
    assert(nodePrcommListWrapper.ptr==NULL);

    // set the pointer and action for consistency checking in Finish
    nodePrcommListWrapper.ptr = (void*) d;
    nodePrcommListWrapper.action = action;

    // do the packing and post the mpi messages...
    updateR2Start(d, action, nodePrcommListWrapper);

  }

  void updateNoR2Finish(double(*d)[3], const int action)
  {

    // consistency check...
    assert(nodePrcommListWrapper.ptr==(void*) d);
    assert(nodePrcommListWrapper.action==action);

    // do the work...
    updateR2Finish(d, action, nodePrcommListWrapper);

    // and reset the pointer...
    nodePrcommListWrapper.ptr = NULL;

  }

  void updateNoR2Reverse(double(*d)[3], const int action)
  {
    updateR2Reverse(d, action, nodePrcommList);
  }

  void updateNoData(double(*d)[3][3], const int action)
  {
    updateR3(d, action, nodePrcommList);
  }

  // (slowly) migrate towards this form...
  void updateNoR3(double(*d)[3][3], const int action)
  {
    updateR3(d, action, nodePrcommList);
  }

  void updateNoSymmetricR3(double(*diag)[3], double(*offd)[3], const int action)
  {
    updateSymmetricR3(diag, offd, action, nodePrcommList);
  }

  // for edges

  void updateEdI1(int * d, const int action)
  {
    updateI1(d, action, edgePrcommList);
  }

  void updateEdI1Reverse(int * d, const int action)
  {
    updateI1Reverse(d, action, edgePrcommList);
  }

  void updateEdI2(int(*d)[3], const int action)
  {
    updateI2(d, action, edgePrcommList);
  }

  void updateEdI2Reverse(int(*d)[3], const int action)
  {
    updateI2Reverse(d, action, edgePrcommList);
  }

  void updateEdData(double * d, const int action)
  {
    updateR1(d, action, edgePrcommList);
  }

  void updateEdR1(double * d, const int action)
  {
    updateR1(d, action, edgePrcommList);
  }

  void updateEdData(double(*d)[3], const int action)
  {
    // make sure the edgePrcommList buffers are not currently
    // in use by another non-blocking message exchange...
    //assert( edgePrcommListWrapper.ptr == NULL ); // something about the non-blocking comm - F. Ham
    updateR2(d, action, edgePrcommList);
  }

  void updateEdR2(double(*d)[3], const int action)
  {
    // make sure the edgePrcommList buffers are not currently
    // in use by another non-blocking message exchange...
    //assert( edgePrcommListWrapper.ptr == NULL ); // something about the non-blocking comm - F. Ham
    updateR2(d, action, edgePrcommList);
  }

  void updateEdR2Reverse(double(*d)[3], const int action)
  {
    updateR2Reverse(d, action, edgePrcommList);
  }

  // general...

protected:

  void updateI1(int * d, const int action, list<Prcomm>& prcommList);

  void updateI1Reverse(int * d, const int action, list<Prcomm>& prcommList);

  void updateI2(int(*d)[3], const int action, list<Prcomm>& prcommList);

  void updateI2Reverse(int(*d)[3], const int action, list<Prcomm>& prcommList);

  void updateR1(double * d, const int action, list<Prcomm>& prcommList);

  void updateR2(double(*d)[3], const int action, list<Prcomm>& prcommList);

  /*
   void updateR1R2R1R1R1R1(double * s1,double (*v1)[3],double * s2,double * s3,double * s4,double * s5,
   const int action,list<Prcomm>& prcommList);
   */

  void updateR3(double(*d)[3][3], const int action, list<Prcomm>& prcommList);

  void updateSymmetricR3(double(*diag)[3], double(*offd)[3], const int action, list<Prcomm>& prcommList);

  void updateSymmetricR3(double(*d)[6], const int action, list<Prcomm>& prcommList);

  void updateR1Start(double *d, const int action, PrcommListWrapper& prcommListWrapper);

  void updateR1Finish(double *d, const int action, PrcommListWrapper& prcommListWrapper);

private:

  void updateR1FinishPrcomm(double *d, const int action, const Prcomm * prcomm);

protected:

  void updateR2Start(double(*d)[3], const int action, PrcommListWrapper& prcommListWrapper);

  void updateR2Finish(double(*d)[3], const int action, PrcommListWrapper& prcommListWrapper);

private:

  void updateR2FinishPrcomm(double(*d)[3], const int action, const Prcomm * prcomm);

protected:

  void updatePackBufferDoublePtrs(list<Prcomm>& prcommList);

  void updateR2Reverse(double(*d)[3], const int action, list<Prcomm>& prcommList);

  void exchangePrcommBufferInt(list<Prcomm>& prcommList)
  {
    exchangePrcommBufferInt(prcommList, 1);
  }

  void exchangePrcommBufferInt(list<Prcomm>& prcommList, const int n);

  void exchangePrcommBufferIntReverse(list<Prcomm>& prcommList)
  {
    exchangePrcommBufferIntReverse(prcommList, 1);
  }

  void exchangePrcommBufferIntReverse(list<Prcomm>& prcommList, const int n);

  void exchangePrcommBufferIntReverseFirst(list<Prcomm>& prcommList, const int n)
  {

    // if this is the first exchange, then only the prcomm->unpackIndexVec is set, and the 
    // packIndexVec should be empty. For certain types, the communicator pairing may not be
    // be symmetric (is this true? - it will error below), so we call this routine to set up 
    // the expected buffer sizes...

    int * send_count = new int[mpi_size];
    for (int i = 0; i<mpi_size; ++i)
      send_count[i] = 0;

    for (list<Prcomm>::iterator prcomm = prcommList.begin(); prcomm!=prcommList.end(); prcomm++)
    {
      assert(prcomm->unpackIndexVec.size()>0);
      assert(prcomm->packIndexVec.size()==0);
      assert(send_count[prcomm->getNbrRank()]==0);
      send_count[prcomm->getNbrRank()] = prcomm->unpackIndexVec.size();
    }

    int * recv_count = new int[mpi_size];
    MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT, mpi_comm);
    delete[] send_count;

    for (int i = 0; i<mpi_size; ++i)
    {
      if (recv_count[i]>0)
      {
        // do we have a Prcomm for this rank?...
        list<Prcomm>::iterator prcomm = prcommList.begin();
        while ((prcomm!=prcommList.end())&&(prcomm->getNbrRank()!=i))
          ++prcomm;
        if (prcomm==prcommList.end())
        {
          // this is an asymetric communicator - it happens someimes
          // with the UgpWithNo node parition. Bascially this is saying the unpack buffer
          // is requesting data from the partition that owns its ghost nodes, but we
          // do not have active nodes that the other processor needs...
          prcommList.push_back(Prcomm(i));
          Prcomm * prcomm = &(prcommList.back());
          assert(prcomm->unpackIndexVec.size()==0);
          assert(prcomm->packIndexVec.size()==0);
          prcomm->packIndexVec.resize(recv_count[i]);
        }
        else
        {
          assert(prcomm->packIndexVec.size()==0);
          prcomm->packIndexVec.resize(recv_count[i]);
        }
      }
    }

    delete[] recv_count;

    // call the normal reverse routine...
    exchangePrcommBufferIntReverse(prcommList, n);

  }

  void exchangePrcommNpackV(list<Prcomm>& prcommList);

  void exchangePrcommNunpackV(list<Prcomm>& prcommList);

  void exchangePrcommBufferIntV(list<Prcomm>& prcommList, const int n);

  void exchangePrcommBufferIntV(list<Prcomm>& prcommList)
  {
    exchangePrcommBufferIntV(prcommList, 1);
  }

  void exchangePrcommBufferDoubleV(list<Prcomm>& prcommList, const int n);

  void exchangePrcommBufferDoubleV(list<Prcomm>& prcommList)
  {
    exchangePrcommBufferDoubleV(prcommList, 1);
  }

  void dumpFacePrcommData();

  void dumpNodePrcommData();

  void dump();

public:

  virtual ~Gp()
  {
  }

};

#endif
