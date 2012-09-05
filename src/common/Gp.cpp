#include "Gp.h"

// mpi tags used for data exchanges...

#define UPDATE_R1_TAG                   1002
#define UPDATE_R2_TAG                   1003
#define UPDATE_R2_REVERSE_TAG           1004
#define EXCHANGE_INT_TAG                1005
#define EXCHANGE_INT_REVERSE_TAG        1006
#define UPDATE_I1_TAG                   1007
#define UPDATE_I1_REVERSE_TAG           1008
#define UPDATE_I2_TAG                   1009
#define UPDATE_I2_REVERSE_TAG           1010
#define UPDATE_R3_TAG                   1011
#define UPDATE_SYMMETRIC_R3_TAG         1012
#define EXCHANGE_PRCOMM_INTV_TAG1       1013
#define EXCHANGE_PRCOMM_INTV_TAG2       1014
#define EXCHANGE_PRCOMM_DOUBLEV_TAG1    1015
#define EXCHANGE_PRCOMM_DOUBLEV_TAG2    1016
#define EXCHANGE_PRCOMM_NPACKV_TAG      1017
#define EXCHANGE_PRCOMM_NUNPACKV_TAG    1018
#define UPDATE_ADDRESS_TAG              1019

// new tags based on RX and IX format...

void Gp::updateR1(double * d, const int action, list<Prcomm>& prcommList)
{

  // count the requests, not including self-requests, which do not 
  // invoke the MPI communication. Requests have to be put into an array
  // ..

  int requestCount = 0;
  list<Prcomm>::iterator prcomm;
  for (prcomm = prcommList.begin(); prcomm!=prcommList.end(); prcomm++)
    if (prcomm->getNbrRank()!=mpi_rank) requestCount++;

  MPI_Request * sendRequestArray;
  MPI_Request * recvRequestArray;
  MPI_Status * statusArray;
  if (requestCount>0)
  {
    sendRequestArray = new MPI_Request[requestCount];
    recvRequestArray = new MPI_Request[requestCount];
    statusArray = new MPI_Status[requestCount];
  }

  // start...

  requestCount = 0;
  for (prcomm = prcommList.begin(); prcomm!=prcommList.end(); prcomm++)
  {

    // ensure size of unpackBufferDouble...

    int nunpack = prcomm->unpackIndexVec.size();
    if (nunpack>prcomm->unpackBufferDoubleSize)
    {
      double * ptr;
      if ((ptr = (double*) realloc(prcomm->unpackBufferDouble, sizeof(double)*nunpack))==NULL)
      {
        cerr<<"Error: realloc failed in updateR1(double * d,..."<<endl;
        throw(-1);
      }
      prcomm->unpackBufferDouble = ptr;
      prcomm->unpackBufferDoubleSize = nunpack;
    }

    // post the Irecv...

    if (prcomm->getNbrRank()!=mpi_rank)
    {
      MPI_Irecv(prcomm->unpackBufferDouble, nunpack, MPI_DOUBLE, prcomm->getNbrRank(), UPDATE_R1_TAG, mpi_comm,
          &(recvRequestArray[requestCount]));
      // XXX we should also store a pointer to the prcomm so when we use
      // MPI_Waitany in the future, we know exactly which prcomm is available
      // for unpacking.
    }

    // ensure size of packBufferDouble...

    int npack = prcomm->packIndexVec.size();
    if (npack>prcomm->packBufferDoubleSize)
    {
      double * ptr;
      if ((ptr = (double*) realloc(prcomm->packBufferDouble, sizeof(double)*npack))==NULL)
      {
        cerr<<"Error: realloc failed in updateR1(double * d,..."<<endl;
        throw(-1);
      }
      prcomm->packBufferDouble = ptr;
      prcomm->packBufferDoubleSize = npack;
    }

    // pack...

    int i;
    for (i = 0; i<npack; i++)
    {
      int index = prcomm->packIndexVec[i];
      prcomm->packBufferDouble[i] = d[index];
    }

    // and send...

    if (prcomm->getNbrRank()==mpi_rank)
    {
      assert(npack==nunpack);
      for (i = 0; i<npack; i++)
        prcomm->unpackBufferDouble[i] = prcomm->packBufferDouble[i];
    }
    else
    {
      MPI_Issend(prcomm->packBufferDouble, npack, MPI_DOUBLE, prcomm->getNbrRank(), UPDATE_R1_TAG, mpi_comm,
          &(sendRequestArray[requestCount]));
      // inc request count...
      requestCount++;
    }
  }

  // could do other work here - so eventually break up this routine into
  // 2 routines - start and finish...

  // now we wait for all messages to be received...

  if (requestCount>0) MPI_Waitall(requestCount, recvRequestArray, statusArray);

  // unpack...

  switch (action)
  {
  case REPLACE_DATA:
    for (prcomm = prcommList.begin(); prcomm!=prcommList.end(); prcomm++)
    {
      int nunpack = prcomm->unpackIndexVec.size();
      int i;
      for (i = 0; i<nunpack; i++)
      {
        int index = prcomm->unpackIndexVec[i];
        d[index] = prcomm->unpackBufferDouble[i];
      }
    }
    break;
  case ADD_DATA:
    for (prcomm = prcommList.begin(); prcomm!=prcommList.end(); prcomm++)
    {
      int nunpack = prcomm->unpackIndexVec.size();
      int i;
      for (i = 0; i<nunpack; i++)
      {
        int index = prcomm->unpackIndexVec[i];
        d[index] += prcomm->unpackBufferDouble[i];
      }
    }
    break;
  case SUBTRACT_DATA:
    for (prcomm = prcommList.begin(); prcomm!=prcommList.end(); prcomm++)
    {
      int nunpack = prcomm->unpackIndexVec.size();
      int i;
      for (i = 0; i<nunpack; i++)
      {
        int index = prcomm->unpackIndexVec[i];
        d[index] -= prcomm->unpackBufferDouble[i];
      }
    }
    break;
  case MIN_DATA:
    for (prcomm = prcommList.begin(); prcomm!=prcommList.end(); prcomm++)
    {
      int nunpack = prcomm->unpackIndexVec.size();
      int i;
      for (i = 0; i<nunpack; i++)
      {
        int index = prcomm->unpackIndexVec[i];
        d[index] = min(d[index], prcomm->unpackBufferDouble[i]);
      }
    }
    break;
  case MAX_DATA:
    for (prcomm = prcommList.begin(); prcomm!=prcommList.end(); prcomm++)
    {
      int nunpack = prcomm->unpackIndexVec.size();
      int i;
      for (i = 0; i<nunpack; i++)
      {
        int index = prcomm->unpackIndexVec[i];
        d[index] = max(d[index], prcomm->unpackBufferDouble[i]);
      }
    }
    break;
  case ADD_NO_PERIODIC_DATA:
    for (prcomm = prcommList.begin(); prcomm!=prcommList.end(); prcomm++)
    {
      list<Range>::iterator range;
      for (range = prcomm->unpackRangeList.begin(); range!=prcomm->unpackRangeList.end(); range++)
      {
        // non-periodic ranges only...
        if (range->isPeriodic()==0)
        {
          for (int i = range->getIndexFirst(); i<=range->getIndexLast(); i++)
          {
            int index = prcomm->unpackIndexVec[i];
            d[index] += prcomm->unpackBufferDouble[i];
          }
        }
      }
    }
    break;
  default:
    cerr<<"Error: unsupported action: "<<action<<endl;
    throw(-1);
  }

  // cleanup after sends...
  if (requestCount>0)
  {
    MPI_Waitall(requestCount, sendRequestArray, statusArray);
    delete[] sendRequestArray;
    delete[] recvRequestArray;
    delete[] statusArray;
  }

}

void Gp::updateR2(double(*d)[3], const int action, list<Prcomm>& prcommList)
{

  // count the requests, not including self-requests (which do not 
  // invoke the MPI communication) Requests have to be put into an array
  // ..

  //MPI_Pause("1");

  int requestCount = 0;
  for (list<Prcomm>::iterator prcomm = prcommList.begin(); prcomm!=prcommList.end(); prcomm++)
    if (prcomm->getNbrRank()!=mpi_rank) requestCount++;

  MPI_Request * sendRequestArray;
  MPI_Request * recvRequestArray;
  MPI_Status * statusArray;
  if (requestCount>0)
  {
    sendRequestArray = new MPI_Request[requestCount];
    recvRequestArray = new MPI_Request[requestCount];
    statusArray = new MPI_Status[requestCount];
  }

  //MPI_Pause("2");

  // start...

  requestCount = 0;
  for (list<Prcomm>::iterator prcomm = prcommList.begin(); prcomm!=prcommList.end(); prcomm++)
  {

    // ensure size of unpackBufferDouble...

    int nunpack = prcomm->unpackIndexVec.size();
    if (3*nunpack>prcomm->unpackBufferDoubleSize)
    {
      double * ptr;
      if ((ptr = (double*) realloc(prcomm->unpackBufferDouble, sizeof(double)*3*nunpack))==NULL)
      {
        cerr<<"Error: realloc failed in updateR2(double (*d)[3],..."<<endl;
        throw(-1);
      }
      prcomm->unpackBufferDouble = ptr;
      prcomm->unpackBufferDoubleSize = 3*nunpack;
    }

    // post the Irecv...

    if (prcomm->getNbrRank()!=mpi_rank)
    {
      MPI_Irecv(prcomm->unpackBufferDouble, 3*nunpack, MPI_DOUBLE, prcomm->getNbrRank(), UPDATE_R2_TAG, mpi_comm,
          &(recvRequestArray[requestCount]));
      // XXX we should also store a pointer to the prcomm so when we use
      // MPI_Waitany in the future, we know exactly which prcomm is available
      // for unpacking.
    }

    //MPI_Pause("3");

    // ensure size of packBufferDouble...

    int npack = prcomm->packIndexVec.size();
    if (3*npack>prcomm->packBufferDoubleSize)
    {
      double * ptr;
      if ((ptr = (double*) realloc(prcomm->packBufferDouble, sizeof(double)*3*npack))==NULL)
      {
        cerr<<"Error: realloc failed in updateR2(double (*d)[3],..."<<endl;
        throw(-1);
      }
      prcomm->packBufferDouble = ptr;
      prcomm->packBufferDoubleSize = 3*npack;
    }

    //MPI_Pause("4");

    // pack...
    switch (action)
    {
    case REPLACE_DATA:
    case ADD_DATA:
    case MIN_DATA:
    case MAX_DATA:
    case ADD_NO_PERIODIC_DATA:
      for (int i = 0; i<npack; i++)
      {
        int index = prcomm->packIndexVec[i];
        for (int j = 0; j<3; j++)
          prcomm->packBufferDouble[3*i+j] = d[index][j];
      }
      break;
    case REPLACE_TRANSLATE_DATA:
    case ADD_TRANSLATE_DATA:
      for (list<Range>::iterator ri = prcomm->packRangeList.begin(); ri!=prcomm->packRangeList.end(); ri++)
      {
        for (int i = ri->getIndexFirst(); i<=ri->getIndexLast(); i++)
        {
          int index = prcomm->packIndexVec[i];
          for (int j = 0; j<3; j++)
            prcomm->packBufferDouble[3*i+j] = d[index][j];
        }
        // rangeTranslateVector is a virtual function that must be defined in 
        // any instantiated version of Gp...
        rangeTranslateVector(prcomm->packBufferDouble+3*ri->getIndexFirst(), // start
            ri->getIndexLast()-ri->getIndexFirst()+1, // size
            ri->getFlag(), ri->dxyz); // info about the transform
      }
      break;
    case REPLACE_ROTATE_DATA:
    case ADD_ROTATE_DATA:
      for (list<Range>::iterator ri = prcomm->packRangeList.begin(); ri!=prcomm->packRangeList.end(); ri++)
      {
        for (int i = ri->getIndexFirst(); i<=ri->getIndexLast(); i++)
        {
          int index = prcomm->packIndexVec[i];
          for (int j = 0; j<3; j++)
            prcomm->packBufferDouble[3*i+j] = d[index][j];
        }
        // rangeRotateVector is a virtual function that must be defined in 
        // any instantiated version of Gp...
        rangeRotateVector(prcomm->packBufferDouble+3*ri->getIndexFirst(), // start
            ri->getIndexLast()-ri->getIndexFirst()+1, // size
            ri->getFlag(), ri->dxyz); // info about the transform
      }
      break;
    default:
      cerr<<"Error: unsupported action: "<<action<<endl;
      throw(-1);
    }

    //MPI_Pause("5");

    // and send...

    if (prcomm->getNbrRank()==mpi_rank)
    {
      assert(npack==nunpack);
      for (int i = 0; i<3*npack; i++)
        prcomm->unpackBufferDouble[i] = prcomm->packBufferDouble[i];
    }
    else
    {
      MPI_Issend(prcomm->packBufferDouble, 3*npack, MPI_DOUBLE, prcomm->getNbrRank(), UPDATE_R2_TAG, mpi_comm,
          &(sendRequestArray[requestCount]));
      // inc request count...
      requestCount++;
    }
  }

  //MPI_Pause("6");


  // could do other work here - so eventually break up this routine into
  // 2 routines - start and finish...

  // now we wait for all messages to be received...

  if (requestCount>0) MPI_Waitall(requestCount, recvRequestArray, statusArray);

  //MPI_Pause("7");


  // unpack...

  switch (action)
  {
  case REPLACE_DATA:
  case REPLACE_TRANSLATE_DATA:
  case REPLACE_ROTATE_DATA:
    for (list<Prcomm>::iterator prcomm = prcommList.begin(); prcomm!=prcommList.end(); prcomm++)
    {
      int nunpack = prcomm->unpackIndexVec.size();
      for (int i = 0; i<nunpack; i++)
      {
        int index = prcomm->unpackIndexVec[i];
        for (int j = 0; j<3; j++)
          d[index][j] = prcomm->unpackBufferDouble[3*i+j];
      }
    }
    break;
  case ADD_DATA:
  case ADD_TRANSLATE_DATA:
  case ADD_ROTATE_DATA:
    for (list<Prcomm>::iterator prcomm = prcommList.begin(); prcomm!=prcommList.end(); prcomm++)
    {
      int nunpack = prcomm->unpackIndexVec.size();
      for (int i = 0; i<nunpack; i++)
      {
        int index = prcomm->unpackIndexVec[i];
        for (int j = 0; j<3; j++)
          d[index][j] += prcomm->unpackBufferDouble[3*i+j];
      }
    }
    break;
  case ADD_NO_PERIODIC_DATA:
    for (list<Prcomm>::iterator prcomm = prcommList.begin(); prcomm!=prcommList.end(); prcomm++)
    {
      list<Range>::iterator range;
      for (range = prcomm->unpackRangeList.begin(); range!=prcomm->unpackRangeList.end(); range++)
      {
        // non-periodic ranges only...
        if (range->isPeriodic()==0)
        {
          for (int i = range->getIndexFirst(); i<=range->getIndexLast(); i++)
          {
            int index = prcomm->unpackIndexVec[i];
            for (int j = 0; j<3; j++)
              d[index][j] += prcomm->unpackBufferDouble[3*i+j];
          }
        }
      }
    }
    break;
  case MIN_DATA:
    for (list<Prcomm>::iterator prcomm = prcommList.begin(); prcomm!=prcommList.end(); prcomm++)
    {
      int nunpack = prcomm->unpackIndexVec.size();
      for (int i = 0; i<nunpack; i++)
      {
        int index = prcomm->unpackIndexVec[i];
        for (int j = 0; j<3; j++)
          d[index][j] = min(d[index][j], prcomm->unpackBufferDouble[3*i+j]);
      }
    }
    break;
  case MAX_DATA:
    for (list<Prcomm>::iterator prcomm = prcommList.begin(); prcomm!=prcommList.end(); prcomm++)
    {
      int nunpack = prcomm->unpackIndexVec.size();
      for (int i = 0; i<nunpack; i++)
      {
        int index = prcomm->unpackIndexVec[i];
        for (int j = 0; j<3; j++)
          d[index][j] = max(d[index][j], prcomm->unpackBufferDouble[3*i+j]);
      }
    }
    break;
  default:
    cerr<<"Error: unsupported action: "<<action<<endl;
    throw(-1);
  }

  // cleanup after sends...
  if (requestCount>0)
  {
    MPI_Waitall(requestCount, sendRequestArray, statusArray);
    delete[] sendRequestArray;
    delete[] recvRequestArray;
    delete[] statusArray;
  }

  //MPI_Pause("8");

}

void Gp::updateR3(double(*d)[3][3], const int action, list<Prcomm>& prcommList)
{

  // count the requests, not including self-requests (which do not 
  // invoke the MPI communication) Requests have to be put into an array
  // ..

  int requestCount = 0;
  for (list<Prcomm>::iterator prcomm = prcommList.begin(); prcomm!=prcommList.end(); prcomm++)
    if (prcomm->getNbrRank()!=mpi_rank) requestCount++;

  MPI_Request * sendRequestArray;
  MPI_Request * recvRequestArray;
  MPI_Status * statusArray;
  if (requestCount>0)
  {
    sendRequestArray = new MPI_Request[requestCount];
    recvRequestArray = new MPI_Request[requestCount];
    statusArray = new MPI_Status[requestCount];
  }

  // start...

  requestCount = 0;
  for (list<Prcomm>::iterator prcomm = prcommList.begin(); prcomm!=prcommList.end(); prcomm++)
  {

    // ensure size of unpackBufferDouble...

    int nunpack = prcomm->unpackIndexVec.size();
    if (9*nunpack>prcomm->unpackBufferDoubleSize)
    {
      double * ptr;
      if ((ptr = (double*) realloc(prcomm->unpackBufferDouble, sizeof(double)*9*nunpack))==NULL)
      {
        cerr<<"Error: realloc failed in updateR3(double (*d)[3][3],..."<<endl;
        throw(-1);
      }
      prcomm->unpackBufferDouble = ptr;
      prcomm->unpackBufferDoubleSize = 9*nunpack;
    }

    // post the Irecv...

    if (prcomm->getNbrRank()!=mpi_rank)
    {
      MPI_Irecv(prcomm->unpackBufferDouble, 9*nunpack, MPI_DOUBLE, prcomm->getNbrRank(), UPDATE_R3_TAG, mpi_comm,
          &(recvRequestArray[requestCount]));
    }

    // ensure size of packBufferDouble...

    int npack = prcomm->packIndexVec.size();
    if (9*npack>prcomm->packBufferDoubleSize)
    {
      double * ptr;
      if ((ptr = (double*) realloc(prcomm->packBufferDouble, sizeof(double)*9*npack))==NULL)
      {
        cerr<<"Error: realloc failed in updateR3(double (*d)[3][3],..."<<endl;
        throw(-1);
      }
      prcomm->packBufferDouble = ptr;
      prcomm->packBufferDoubleSize = 9*npack;
    }

    // pack...
    switch (action)
    {
    case REPLACE_ROTATE_DATA:
    case ADD_ROTATE_DATA:
    case SUBTRACT_ROTATE_DATA:
      for (list<Range>::iterator ri = prcomm->packRangeList.begin(); ri!=prcomm->packRangeList.end(); ri++)
      {
        for (int i = ri->getIndexFirst(); i<=ri->getIndexLast(); i++)
        {
          int index = prcomm->packIndexVec[i];
          for (int j = 0; j<3; j++)
            for (int k = 0; k<3; k++)
              prcomm->packBufferDouble[9*i+3*j+k] = d[index][j][k];
        }
        // rangeRotateTensor is a virtual function that must be defined in 
        // any instantiated version of Gp...
        rangeRotateTensor(prcomm->packBufferDouble+9*ri->getIndexFirst(), // start
            ri->getIndexLast()-ri->getIndexFirst()+1, // size
            ri->getFlag(), ri->dxyz); // info about the transform
      }
      break;
    default:
      cerr<<"Error: unsupported action: "<<action<<endl;
      throw(-1);
    }

    // and send...

    if (prcomm->getNbrRank()==mpi_rank)
    {
      assert(npack==nunpack);
      for (int i = 0; i<9*npack; i++)
        prcomm->unpackBufferDouble[i] = prcomm->packBufferDouble[i];
    }
    else
    {
      MPI_Issend(prcomm->packBufferDouble, 9*npack, MPI_DOUBLE, prcomm->getNbrRank(), UPDATE_R3_TAG, mpi_comm,
          &(sendRequestArray[requestCount]));
      // inc request count...
      requestCount++;
    }
  }

  // could do other work here - so eventually break up this routine into
  // 2 routines - start and finish...

  // now we wait for all messages to be received...

  if (requestCount>0) MPI_Waitall(requestCount, recvRequestArray, statusArray);

  // unpack...

  switch (action)
  {
  case REPLACE_ROTATE_DATA:
    for (list<Prcomm>::iterator prcomm = prcommList.begin(); prcomm!=prcommList.end(); prcomm++)
    {
      int nunpack = prcomm->unpackIndexVec.size();
      for (int i = 0; i<nunpack; i++)
      {
        int index = prcomm->unpackIndexVec[i];
        for (int j = 0; j<3; j++)
          for (int k = 0; k<3; k++)
            d[index][j][k] = prcomm->unpackBufferDouble[9*i+3*j+k];
      }
    }
    break;
  case ADD_ROTATE_DATA:
    for (list<Prcomm>::iterator prcomm = prcommList.begin(); prcomm!=prcommList.end(); prcomm++)
    {
      int nunpack = prcomm->unpackIndexVec.size();
      for (int i = 0; i<nunpack; i++)
      {
        int index = prcomm->unpackIndexVec[i];
        for (int j = 0; j<3; j++)
          for (int k = 0; k<3; k++)
            d[index][j][k] += prcomm->unpackBufferDouble[9*i+3*j+k];
      }
    }
    break;
  case SUBTRACT_ROTATE_DATA:
    for (list<Prcomm>::iterator prcomm = prcommList.begin(); prcomm!=prcommList.end(); prcomm++)
    {
      int nunpack = prcomm->unpackIndexVec.size();
      for (int i = 0; i<nunpack; i++)
      {
        int index = prcomm->unpackIndexVec[i];
        for (int j = 0; j<3; j++)
          for (int k = 0; k<3; k++)
            d[index][j][k] -= prcomm->unpackBufferDouble[9*i+3*j+k];
      }
    }
    break;
  default:
    cerr<<"Error: unsupported action: "<<action<<endl;
    throw(-1);
  }

  // cleanup after sends...
  if (requestCount>0)
  {
    MPI_Waitall(requestCount, sendRequestArray, statusArray);
    delete[] sendRequestArray;
    delete[] recvRequestArray;
    delete[] statusArray;
  }

}

void Gp::updateSymmetricR3(double(*diag)[3], double(*offd)[3], const int action, list<Prcomm>& prcommList)
{

  // count the requests, not including self-requests (which do not 
  // invoke the MPI communication) Requests have to be put into an array
  // ..

  int requestCount = 0;
  for (list<Prcomm>::iterator prcomm = prcommList.begin(); prcomm!=prcommList.end(); prcomm++)
    if (prcomm->getNbrRank()!=mpi_rank) requestCount++;

  MPI_Request * sendRequestArray;
  MPI_Request * recvRequestArray;
  MPI_Status * statusArray;
  if (requestCount>0)
  {
    sendRequestArray = new MPI_Request[requestCount];
    recvRequestArray = new MPI_Request[requestCount];
    statusArray = new MPI_Status[requestCount];
  }

  // start...

  requestCount = 0;
  for (list<Prcomm>::iterator prcomm = prcommList.begin(); prcomm!=prcommList.end(); prcomm++)
  {

    // ensure size of unpackBufferDouble...

    int nunpack = prcomm->unpackIndexVec.size();
    if (6*nunpack>prcomm->unpackBufferDoubleSize)
    {
      double * ptr;
      if ((ptr = (double*) realloc(prcomm->unpackBufferDouble, sizeof(double)*6*nunpack))==NULL)
      {
        cerr<<"Error: realloc failed in updateSymmetricR3(double (*diag)[3],..."<<endl;
        throw(-1);
      }
      prcomm->unpackBufferDouble = ptr;
      prcomm->unpackBufferDoubleSize = 6*nunpack;
    }

    // post the Irecv...

    if (prcomm->getNbrRank()!=mpi_rank)
    {
      MPI_Irecv(prcomm->unpackBufferDouble, 6*nunpack, MPI_DOUBLE, prcomm->getNbrRank(), UPDATE_SYMMETRIC_R3_TAG,
          mpi_comm, &(recvRequestArray[requestCount]));
    }

    // ensure size of packBufferDouble...

    int npack = prcomm->packIndexVec.size();
    if (6*npack>prcomm->packBufferDoubleSize)
    {
      double * ptr;
      if ((ptr = (double*) realloc(prcomm->packBufferDouble, sizeof(double)*6*npack))==NULL)
      {
        cerr<<"Error: realloc failed in updateSymmetricR3(double (*diag)[3],..."<<endl;
        throw(-1);
      }
      prcomm->packBufferDouble = ptr;
      prcomm->packBufferDoubleSize = 6*npack;
    }

    // pack...
    switch (action)
    {
    case REPLACE_ROTATE_DATA:
    case ADD_ROTATE_DATA:
      for (list<Range>::iterator ri = prcomm->packRangeList.begin(); ri!=prcomm->packRangeList.end(); ri++)
      {
        for (int i = ri->getIndexFirst(); i<=ri->getIndexLast(); i++)
        {
          int index = prcomm->packIndexVec[i];
          for (int j = 0; j<3; j++)
            prcomm->packBufferDouble[6*i+j] = diag[index][j];
          for (int j = 0; j<3; j++)
            prcomm->packBufferDouble[6*i+3+j] = offd[index][j];
        }
        // rangeRotateSymmetricTensor is a virtual function that must be defined in 
        // any instantiated version of Gp...
        rangeRotateSymmetricTensor(prcomm->packBufferDouble+6*ri->getIndexFirst(), // start
            ri->getIndexLast()-ri->getIndexFirst()+1, // size
            ri->getFlag(), ri->dxyz); // info about the transform
      }
      break;
    default:
      cerr<<"Error: unsupported action: "<<action<<endl;
      throw(-1);
    }

    // and send...

    if (prcomm->getNbrRank()==mpi_rank)
    {
      assert(npack==nunpack);
      for (int i = 0; i<6*npack; i++)
        prcomm->unpackBufferDouble[i] = prcomm->packBufferDouble[i];
    }
    else
    {
      MPI_Issend(prcomm->packBufferDouble, 6*npack, MPI_DOUBLE, prcomm->getNbrRank(), UPDATE_SYMMETRIC_R3_TAG,
          mpi_comm, &(sendRequestArray[requestCount]));
      // inc request count...
      requestCount++;
    }
  }

  // could do other work here - so eventually break up this routine into
  // 2 routines - start and finish...

  // now we wait for all messages to be received...

  if (requestCount>0) MPI_Waitall(requestCount, recvRequestArray, statusArray);

  // unpack...

  switch (action)
  {
  case REPLACE_ROTATE_DATA:
    for (list<Prcomm>::iterator prcomm = prcommList.begin(); prcomm!=prcommList.end(); prcomm++)
    {
      int nunpack = prcomm->unpackIndexVec.size();
      for (int i = 0; i<nunpack; i++)
      {
        int index = prcomm->unpackIndexVec[i];
        for (int j = 0; j<3; j++)
          diag[index][j] = prcomm->unpackBufferDouble[6*i+j];
        for (int j = 0; j<3; j++)
          offd[index][j] = prcomm->unpackBufferDouble[6*i+3+j];
      }
    }
    break;
  case ADD_ROTATE_DATA:
    for (list<Prcomm>::iterator prcomm = prcommList.begin(); prcomm!=prcommList.end(); prcomm++)
    {
      int nunpack = prcomm->unpackIndexVec.size();
      for (int i = 0; i<nunpack; i++)
      {
        int index = prcomm->unpackIndexVec[i];
        for (int j = 0; j<3; j++)
          diag[index][j] += prcomm->unpackBufferDouble[6*i+j];
        for (int j = 0; j<3; j++)
          offd[index][j] += prcomm->unpackBufferDouble[6*i+3+j];
      }
    }
    break;
  default:
    cerr<<"Error: unsupported action: "<<action<<endl;
    throw(-1);
  }

  // cleanup after sends...
  if (requestCount>0)
  {
    MPI_Waitall(requestCount, sendRequestArray, statusArray);
    delete[] sendRequestArray;
    delete[] recvRequestArray;
    delete[] statusArray;
  }

}

void Gp::updateSymmetricR3(double(*d)[6], const int action, list<Prcomm>& prcommList)
{

  // count the requests, not including self-requests (which do not 
  // invoke the MPI communication) Requests have to be put into an array
  // ..

  int requestCount = 0;
  for (list<Prcomm>::iterator prcomm = prcommList.begin(); prcomm!=prcommList.end(); prcomm++)
    if (prcomm->getNbrRank()!=mpi_rank) requestCount++;

  MPI_Request * sendRequestArray;
  MPI_Request * recvRequestArray;
  MPI_Status * statusArray;
  if (requestCount>0)
  {
    sendRequestArray = new MPI_Request[requestCount];
    recvRequestArray = new MPI_Request[requestCount];
    statusArray = new MPI_Status[requestCount];
  }

  // start...

  requestCount = 0;
  for (list<Prcomm>::iterator prcomm = prcommList.begin(); prcomm!=prcommList.end(); prcomm++)
  {

    // ensure size of unpackBufferDouble...

    int nunpack = prcomm->unpackIndexVec.size();
    if (6*nunpack>prcomm->unpackBufferDoubleSize)
    {
      double * ptr;
      if ((ptr = (double*) realloc(prcomm->unpackBufferDouble, sizeof(double)*6*nunpack))==NULL)
      {
        cerr<<"Error: realloc failed in updateSymmetricR3(double (*d)[6],..."<<endl;
        throw(-1);
      }
      prcomm->unpackBufferDouble = ptr;
      prcomm->unpackBufferDoubleSize = 6*nunpack;
    }

    // post the Irecv...

    if (prcomm->getNbrRank()!=mpi_rank)
    {
      MPI_Irecv(prcomm->unpackBufferDouble, 6*nunpack, MPI_DOUBLE, prcomm->getNbrRank(), UPDATE_SYMMETRIC_R3_TAG,
          mpi_comm, &(recvRequestArray[requestCount]));
    }

    // ensure size of packBufferDouble...

    int npack = prcomm->packIndexVec.size();
    if (6*npack>prcomm->packBufferDoubleSize)
    {
      double * ptr;
      if ((ptr = (double*) realloc(prcomm->packBufferDouble, sizeof(double)*6*npack))==NULL)
      {
        cerr<<"Error: realloc failed in updateSymmetricR3(double (*d)[6],..."<<endl;
        throw(-1);
      }
      prcomm->packBufferDouble = ptr;
      prcomm->packBufferDoubleSize = 6*npack;
    }

    // pack...
    switch (action)
    {
    case REPLACE_ROTATE_DATA:
    case ADD_ROTATE_DATA:
      for (list<Range>::iterator ri = prcomm->packRangeList.begin(); ri!=prcomm->packRangeList.end(); ri++)
      {
        for (int i = ri->getIndexFirst(); i<=ri->getIndexLast(); i++)
        {
          int index = prcomm->packIndexVec[i];
          // recall the 6 tensor elements are stored as xx,xy,xz,yy,yz,zz
          prcomm->packBufferDouble[6*i] = d[index][0]; // xx
          prcomm->packBufferDouble[6*i+1] = d[index][3]; // yy
          prcomm->packBufferDouble[6*i+2] = d[index][5]; // zz
          prcomm->packBufferDouble[6*i+3] = d[index][4]; // yz
          prcomm->packBufferDouble[6*i+4] = d[index][2]; // xz
          prcomm->packBufferDouble[6*i+5] = d[index][1]; // xy
        }
        // rangeRotateSymmetricTensor is a virtual function that must be defined in 
        // any instantiated version of Gp...
        rangeRotateSymmetricTensor(prcomm->packBufferDouble+6*ri->getIndexFirst(), // start
            ri->getIndexLast()-ri->getIndexFirst()+1, // size
            ri->getFlag(), ri->dxyz); // info about the transform
      }
      break;
    default:
      cerr<<"Error: unsupported action: "<<action<<endl;
      throw(-1);
    }

    // and send...

    if (prcomm->getNbrRank()==mpi_rank)
    {
      assert(npack==nunpack);
      for (int i = 0; i<6*npack; i++)
        prcomm->unpackBufferDouble[i] = prcomm->packBufferDouble[i];
    }
    else
    {
      MPI_Issend(prcomm->packBufferDouble, 6*npack, MPI_DOUBLE, prcomm->getNbrRank(), UPDATE_SYMMETRIC_R3_TAG,
          mpi_comm, &(sendRequestArray[requestCount]));
      // inc request count...
      requestCount++;
    }
  }

  // could do other work here - so eventually break up this routine into
  // 2 routines - start and finish...

  // now we wait for all messages to be received...

  if (requestCount>0) MPI_Waitall(requestCount, recvRequestArray, statusArray);

  // unpack...

  switch (action)
  {
  case REPLACE_ROTATE_DATA:
    for (list<Prcomm>::iterator prcomm = prcommList.begin(); prcomm!=prcommList.end(); prcomm++)
    {
      int nunpack = prcomm->unpackIndexVec.size();
      for (int i = 0; i<nunpack; i++)
      {
        int index = prcomm->unpackIndexVec[i];
        d[index][0] = prcomm->unpackBufferDouble[6*i]; // xx
        d[index][1] = prcomm->unpackBufferDouble[6*i+5]; // xy
        d[index][2] = prcomm->unpackBufferDouble[6*i+4]; // xz
        d[index][3] = prcomm->unpackBufferDouble[6*i+1]; // yy
        d[index][4] = prcomm->unpackBufferDouble[6*i+3]; // yz
        d[index][5] = prcomm->unpackBufferDouble[6*i+2]; // zz
      }
    }
    break;
  case ADD_ROTATE_DATA:
    for (list<Prcomm>::iterator prcomm = prcommList.begin(); prcomm!=prcommList.end(); prcomm++)
    {
      int nunpack = prcomm->unpackIndexVec.size();
      for (int i = 0; i<nunpack; i++)
      {
        int index = prcomm->unpackIndexVec[i];
        d[index][0] += prcomm->unpackBufferDouble[6*i]; // xx
        d[index][1] += prcomm->unpackBufferDouble[6*i+5]; // xy
        d[index][2] += prcomm->unpackBufferDouble[6*i+4]; // xz
        d[index][3] += prcomm->unpackBufferDouble[6*i+1]; // yy
        d[index][4] += prcomm->unpackBufferDouble[6*i+3]; // yz
        d[index][5] += prcomm->unpackBufferDouble[6*i+2]; // zz
      }
    }
    break;
  default:
    cerr<<"Error: unsupported action: "<<action<<endl;
    throw(-1);
  }

  // cleanup after sends...
  if (requestCount>0)
  {
    MPI_Waitall(requestCount, sendRequestArray, statusArray);
    delete[] sendRequestArray;
    delete[] recvRequestArray;
    delete[] statusArray;
  }

}

void Gp::updateR1Start(double *d, const int action, PrcommListWrapper& prcommListWrapper)
{

  // XXXXXX this is really bad form, but just trying to fix...

  /*
   static int first = 1;
   if (first) {
   cout << "updatePackBufferDoublePtrs" << endl;
   first = 0;
   updatePackBufferDoublePtrs(*prcommListWrapper.prcommListPtr);
   }
   */

  // initialize the wrapper if it has not been initialized yet...

  if (prcommListWrapper.nonSelfPeriodicSize==-1)
  {
    prcommListWrapper.nonSelfPeriodicSize = 0;
    for (list<Prcomm>::iterator prcomm = prcommListWrapper.prcommListPtr->begin(); prcomm
        !=prcommListWrapper.prcommListPtr->end(); prcomm++)
      if (prcomm->getNbrRank()!=mpi_rank) prcommListWrapper.nonSelfPeriodicSize += 1;
    prcommListWrapper.sendRequestArray = new MPI_Request[prcommListWrapper.nonSelfPeriodicSize];
    prcommListWrapper.recvRequestArray = new MPI_Request[prcommListWrapper.nonSelfPeriodicSize];
    prcommListWrapper.statusArray = new MPI_Status[prcommListWrapper.nonSelfPeriodicSize];
    // and set a back-ptr to the appropriate prcomm member in the list...
    prcommListWrapper.prcommPtrVec.resize(prcommListWrapper.nonSelfPeriodicSize);
    int requestCount = 0;
    for (list<Prcomm>::iterator prcomm = prcommListWrapper.prcommListPtr->begin(); prcomm
        !=prcommListWrapper.prcommListPtr->end(); prcomm++)
      if (prcomm->getNbrRank()!=mpi_rank) prcommListWrapper.prcommPtrVec[requestCount++] = &(*prcomm);
  }

  {
    // start...
    int requestCount = 0;
    for (list<Prcomm>::iterator prcomm = prcommListWrapper.prcommListPtr->begin(); prcomm
        !=prcommListWrapper.prcommListPtr->end(); prcomm++)
    {

      // ensure size of unpackBufferDouble...

      int nunpack = prcomm->unpackIndexVec.size();
      if (nunpack>prcomm->unpackBufferDoubleSize)
      {
        double * ptr;
        if ((ptr = (double*) realloc(prcomm->unpackBufferDouble, sizeof(double)*nunpack))==NULL)
        {
          cerr<<"Error: realloc failed in updateR1Start(double *d,..."<<endl;
          throw(-1);
        }
        prcomm->unpackBufferDouble = ptr;
        prcomm->unpackBufferDoubleSize = nunpack;
      }

      // post the Irecv...

      if (prcomm->getNbrRank()!=mpi_rank)
      {
        MPI_Irecv(prcomm->unpackBufferDouble, nunpack, MPI_DOUBLE, prcomm->getNbrRank(), UPDATE_R1_TAG, mpi_comm,
            &(prcommListWrapper.recvRequestArray[requestCount]));
      }

      // ensure size of packBufferDouble...

      int npack = prcomm->packIndexVec.size();
      if (npack>prcomm->packBufferDoubleSize)
      {
        double * ptr;
        if ((ptr = (double*) realloc(prcomm->packBufferDouble, sizeof(double)*npack))==NULL)
        {
          cerr<<"Error: realloc failed in updateR1Start(double *d,..."<<endl;
          throw(-1);
        }
        prcomm->packBufferDouble = ptr;
        prcomm->packBufferDoubleSize = npack;
      }

      /*
       cout << "mpi_rank: " << mpi_rank << " prcomm->getNbrRank(): " << prcomm->getNbrRank() <<
       " prcomm->packBufferDouble: " << prcomm->packBufferDouble << endl;
       MPI_Pause("XXXXX");
       */

      /*
       cout << "mpi_rank: " << mpi_rank << " prcomm->getNbrRank(): " << prcomm->getNbrRank() <<
       " prcomm->packBufferDouble: " << prcomm->packBufferDouble << endl;
       */

      // pack...
      switch (action)
      {
      case REPLACE_DATA:
      case ADD_DATA:
      case MIN_DATA:
      case MAX_DATA:
        for (int i = 0; i<npack; i++)
        {
          int index = prcomm->packIndexVec[i];
          prcomm->packBufferDouble[i] = d[index];
        }
        break;
      default:
        cerr<<"Error: unsupported action: "<<action<<endl;
        throw(-1);
      }

      // and post the send...

      if (prcomm->getNbrRank()!=mpi_rank)
      {
        MPI_Issend(prcomm->packBufferDouble, npack, MPI_DOUBLE, prcomm->getNbrRank(), UPDATE_R1_TAG, mpi_comm,
            &(prcommListWrapper.sendRequestArray[requestCount]));
        // inc request count...
        requestCount++;
      }
    }
  }

  //MPI_Pause("VVV");

}

void Gp::updateR2Start(double(*d)[3], const int action, PrcommListWrapper& prcommListWrapper)
{

  // XXXXXX this is really bad form, but just trying to fix...

  /*
   static int first = 1;
   if (first) {
   cout << "updatePackBufferDoublePtrs" << endl;
   first = 0;
   updatePackBufferDoublePtrs(*prcommListWrapper.prcommListPtr);
   }
   */

  // initialize the wrapper if it has not been initialized yet...

  if (prcommListWrapper.nonSelfPeriodicSize==-1)
  {
    prcommListWrapper.nonSelfPeriodicSize = 0;
    for (list<Prcomm>::iterator prcomm = prcommListWrapper.prcommListPtr->begin(); prcomm
        !=prcommListWrapper.prcommListPtr->end(); prcomm++)
      if (prcomm->getNbrRank()!=mpi_rank) prcommListWrapper.nonSelfPeriodicSize += 1;
    prcommListWrapper.sendRequestArray = new MPI_Request[prcommListWrapper.nonSelfPeriodicSize];
    prcommListWrapper.recvRequestArray = new MPI_Request[prcommListWrapper.nonSelfPeriodicSize];
    prcommListWrapper.statusArray = new MPI_Status[prcommListWrapper.nonSelfPeriodicSize];
    // and set a back-ptr to the appropriate prcomm member in the list...
    prcommListWrapper.prcommPtrVec.resize(prcommListWrapper.nonSelfPeriodicSize);
    int requestCount = 0;
    for (list<Prcomm>::iterator prcomm = prcommListWrapper.prcommListPtr->begin(); prcomm
        !=prcommListWrapper.prcommListPtr->end(); prcomm++)
      if (prcomm->getNbrRank()!=mpi_rank) prcommListWrapper.prcommPtrVec[requestCount++] = &(*prcomm);
  }

  {
    // start...
    int requestCount = 0;
    for (list<Prcomm>::iterator prcomm = prcommListWrapper.prcommListPtr->begin(); prcomm
        !=prcommListWrapper.prcommListPtr->end(); prcomm++)
    {

      // ensure size of unpackBufferDouble...

      int nunpack = prcomm->unpackIndexVec.size();
      if (3*nunpack>prcomm->unpackBufferDoubleSize)
      {
        double * ptr;
        if ((ptr = (double*) realloc(prcomm->unpackBufferDouble, sizeof(double)*3*nunpack))==NULL)
        {
          cerr<<"Error: realloc failed in updateR2Start(double (*d)[3],..."<<endl;
          throw(-1);
        }
        prcomm->unpackBufferDouble = ptr;
        prcomm->unpackBufferDoubleSize = 3*nunpack;
      }

      // post the Irecv...

      if (prcomm->getNbrRank()!=mpi_rank)
      {
        MPI_Irecv(prcomm->unpackBufferDouble, 3*nunpack, MPI_DOUBLE, prcomm->getNbrRank(), UPDATE_R2_TAG, mpi_comm,
            &(prcommListWrapper.recvRequestArray[requestCount]));
      }

      // ensure size of packBufferDouble...

      int npack = prcomm->packIndexVec.size();
      if (3*npack>prcomm->packBufferDoubleSize)
      {
        double * ptr;
        if ((ptr = (double*) realloc(prcomm->packBufferDouble, sizeof(double)*3*npack))==NULL)
        {
          cerr<<"Error: realloc failed in updateR2Start(double (*d)[3],..."<<endl;
          throw(-1);
        }
        prcomm->packBufferDouble = ptr;
        prcomm->packBufferDoubleSize = 3*npack;
      }

      /*
       cout << "mpi_rank: " << mpi_rank << " prcomm->getNbrRank(): " << prcomm->getNbrRank() <<
       " prcomm->packBufferDouble: " << prcomm->packBufferDouble << endl;
       MPI_Pause("XXXXX");
       */

      /*
       cout << "mpi_rank: " << mpi_rank << " prcomm->getNbrRank(): " << prcomm->getNbrRank() <<
       " prcomm->packBufferDouble: " << prcomm->packBufferDouble << endl;
       */

      // pack...
      switch (action)
      {
      case REPLACE_DATA:
      case ADD_DATA:
      case MIN_DATA:
      case MAX_DATA:
        for (int i = 0; i<npack; i++)
        {
          int index = prcomm->packIndexVec[i];
          for (int j = 0; j<3; j++)
            prcomm->packBufferDouble[3*i+j] = d[index][j];
        }
        break;
      case REPLACE_TRANSLATE_DATA:
      case ADD_TRANSLATE_DATA:
        for (list<Range>::iterator ri = prcomm->packRangeList.begin(); ri!=prcomm->packRangeList.end(); ri++)
        {
          for (int i = ri->getIndexFirst(); i<=ri->getIndexLast(); i++)
          {
            int index = prcomm->packIndexVec[i];
            for (int j = 0; j<3; j++)
              prcomm->packBufferDouble[3*i+j] = d[index][j];
          }
          // rangeTranslateVector is a virtual function that must be defined in
          // any instantiated version of Gp...
          rangeTranslateVector(prcomm->packBufferDouble+3*ri->getIndexFirst(), // start
              ri->getIndexLast()-ri->getIndexFirst()+1, // size
              ri->getFlag(), ri->dxyz); // info about the transform
        }
        break;
      case REPLACE_ROTATE_DATA:
      case ADD_ROTATE_DATA:
        for (list<Range>::iterator ri = prcomm->packRangeList.begin(); ri!=prcomm->packRangeList.end(); ri++)
        {
          for (int i = ri->getIndexFirst(); i<=ri->getIndexLast(); i++)
          {
            int index = prcomm->packIndexVec[i];
            for (int j = 0; j<3; j++)
              prcomm->packBufferDouble[3*i+j] = d[index][j];
          }
          // rangeRotateVector is a virtual function that must be defined in
          // any instantiated version of Gp...
          rangeRotateVector(prcomm->packBufferDouble+3*ri->getIndexFirst(), // start
              ri->getIndexLast()-ri->getIndexFirst()+1, // size
              ri->getFlag(), ri->dxyz); // info about the transform
        }
        break;
      default:
        cerr<<"Error: unsupported action: "<<action<<endl;
        throw(-1);
      }

      // and post the send...

      if (prcomm->getNbrRank()!=mpi_rank)
      {
        MPI_Issend(prcomm->packBufferDouble, 3*npack, MPI_DOUBLE, prcomm->getNbrRank(), UPDATE_R2_TAG, mpi_comm,
            &(prcommListWrapper.sendRequestArray[requestCount]));
        // inc request count...
        requestCount++;
      }
    }
  }

  //MPI_Pause("VVV");

}

void Gp::updateR1Finish(double *d, const int action, PrcommListWrapper& prcommListWrapper)
{

  //MPI_Barrier(mpi_comm);

  // self-periodic part first...

  for (list<Prcomm>::iterator prcomm = prcommListWrapper.prcommListPtr->begin(); prcomm
      !=prcommListWrapper.prcommListPtr->end(); prcomm++)
  {
    if (prcomm->getNbrRank()==mpi_rank)
    {
      int nunpack = prcomm->unpackIndexVec.size();
      int npack = prcomm->packIndexVec.size();
      assert(npack==nunpack);
      for (int i = 0; i<npack; i++)
        prcomm->unpackBufferDouble[i] = prcomm->packBufferDouble[i];
      updateR1FinishPrcomm(d, action, &(*prcomm));
    }
  }

  // complete non-blocking mpi part...

  int requestCount = prcommListWrapper.nonSelfPeriodicSize;
  while (requestCount>0)
  {
    int index;
    MPI_Status status;
    MPI_Waitany(prcommListWrapper.nonSelfPeriodicSize, prcommListWrapper.recvRequestArray, &index, &status);
    assert((index>=0)&&(index<prcommListWrapper.nonSelfPeriodicSize));
    updateR1FinishPrcomm(d, action, prcommListWrapper.prcommPtrVec[index]);
    requestCount--;
  }

  // cleanup sends...
  if (prcommListWrapper.nonSelfPeriodicSize>0) MPI_Waitall(prcommListWrapper.nonSelfPeriodicSize,
      prcommListWrapper.sendRequestArray, prcommListWrapper.statusArray);

  //MPI_Pause("WWW");

}

void Gp::updateR1FinishPrcomm(double *d, const int action, const Prcomm * prcomm)
{

  /*
   cout << "mpi_rank: " << mpi_rank << " prcomm->getNbrRank(): " << prcomm->getNbrRank() <<
   " prcomm->packBufferDoublePtr: " << prcomm->packBufferDoublePtr << endl;
   */

  /*
   for (int i = 0; i < min((size_t)4,prcomm->unpackIndexVec.size()); i++) {
   for (int j = 0; j < 3; j++)
   cout << "got unpack, shmem: " << prcomm->unpackBufferDouble[3*i+j] << " " <<
   prcomm->packBufferDoublePtr[3*i+j] << endl;
   }
   */

  switch (action)
  {
  case REPLACE_DATA:
  {
    int nunpack = prcomm->unpackIndexVec.size();
    for (int i = 0; i<nunpack; i++)
    {
      int index = prcomm->unpackIndexVec[i];
      d[index] = prcomm->unpackBufferDouble[i];
    }
  }
    break;
  case ADD_DATA:
  {
    int nunpack = prcomm->unpackIndexVec.size();
    for (int i = 0; i<nunpack; i++)
    {
      int index = prcomm->unpackIndexVec[i];
      d[index] += prcomm->unpackBufferDouble[i];
    }
  }
    break;
  case MIN_DATA:
  {
    int nunpack = prcomm->unpackIndexVec.size();
    for (int i = 0; i<nunpack; i++)
    {
      int index = prcomm->unpackIndexVec[i];
      d[index] = min(d[index], prcomm->unpackBufferDouble[i]);
    }
  }
    break;
  case MAX_DATA:
  {
    int nunpack = prcomm->unpackIndexVec.size();
    for (int i = 0; i<nunpack; i++)
    {
      int index = prcomm->unpackIndexVec[i];
      d[index] = max(d[index], prcomm->unpackBufferDouble[i]);
    }
  }
    break;
  default:
    cerr<<"Error: unsupported action: "<<action<<endl;
    throw(-1);
  }

}

void Gp::updateR2Finish(double(*d)[3], const int action, PrcommListWrapper& prcommListWrapper)
{

  //MPI_Barrier(mpi_comm);

  // self-periodic part first...

  for (list<Prcomm>::iterator prcomm = prcommListWrapper.prcommListPtr->begin(); prcomm
      !=prcommListWrapper.prcommListPtr->end(); prcomm++)
  {
    if (prcomm->getNbrRank()==mpi_rank)
    {
      int nunpack = prcomm->unpackIndexVec.size();
      int npack = prcomm->packIndexVec.size();
      assert(npack==nunpack);
      for (int i = 0; i<3*npack; i++)
        prcomm->unpackBufferDouble[i] = prcomm->packBufferDouble[i];
      updateR2FinishPrcomm(d, action, &(*prcomm));
    }
  }

  // complete non-blocking mpi part...

  int requestCount = prcommListWrapper.nonSelfPeriodicSize;
  while (requestCount>0)
  {
    int index;
    MPI_Status status;
    MPI_Waitany(prcommListWrapper.nonSelfPeriodicSize, prcommListWrapper.recvRequestArray, &index, &status);
    assert((index>=0)&&(index<prcommListWrapper.nonSelfPeriodicSize));
    updateR2FinishPrcomm(d, action, prcommListWrapper.prcommPtrVec[index]);
    requestCount--;
  }

  // cleanup sends...
  if (prcommListWrapper.nonSelfPeriodicSize>0) MPI_Waitall(prcommListWrapper.nonSelfPeriodicSize,
      prcommListWrapper.sendRequestArray, prcommListWrapper.statusArray);

  //MPI_Pause("WWW");

}

void Gp::updateR2FinishPrcomm(double(*d)[3], const int action, const Prcomm * prcomm)
{

  /*
   cout << "mpi_rank: " << mpi_rank << " prcomm->getNbrRank(): " << prcomm->getNbrRank() <<
   " prcomm->packBufferDoublePtr: " << prcomm->packBufferDoublePtr << endl;
   */

  /*
   for (int i = 0; i < min((size_t)4,prcomm->unpackIndexVec.size()); i++) {
   for (int j = 0; j < 3; j++)
   cout << "got unpack, shmem: " << prcomm->unpackBufferDouble[3*i+j] << " " <<
   prcomm->packBufferDoublePtr[3*i+j] << endl;
   }
   */

  switch (action)
  {
  case REPLACE_DATA:
  case REPLACE_TRANSLATE_DATA:
  case REPLACE_ROTATE_DATA:
  {
    int nunpack = prcomm->unpackIndexVec.size();
    for (int i = 0; i<nunpack; i++)
    {
      int index = prcomm->unpackIndexVec[i];
      for (int j = 0; j<3; j++)
        d[index][j] = prcomm->unpackBufferDouble[3*i+j];
    }
  }
    break;
  case ADD_DATA:
  case ADD_TRANSLATE_DATA:
  case ADD_ROTATE_DATA:
  {
    int nunpack = prcomm->unpackIndexVec.size();
    for (int i = 0; i<nunpack; i++)
    {
      int index = prcomm->unpackIndexVec[i];
      for (int j = 0; j<3; j++)
        d[index][j] += prcomm->unpackBufferDouble[3*i+j];
    }
  }
    break;
  case MIN_DATA:
  {
    int nunpack = prcomm->unpackIndexVec.size();
    for (int i = 0; i<nunpack; i++)
    {
      int index = prcomm->unpackIndexVec[i];
      for (int j = 0; j<3; j++)
        d[index][j] = min(d[index][j], prcomm->unpackBufferDouble[3*i+j]);
    }
  }
    break;
  case MAX_DATA:
  {
    int nunpack = prcomm->unpackIndexVec.size();
    for (int i = 0; i<nunpack; i++)
    {
      int index = prcomm->unpackIndexVec[i];
      for (int j = 0; j<3; j++)
        d[index][j] = max(d[index][j], prcomm->unpackBufferDouble[3*i+j]);
    }
  }
    break;
  default:
    cerr<<"Error: unsupported action: "<<action<<endl;
    throw(-1);
  }

}

void Gp::updateR2Reverse(double(*d)[3], const int action, list<Prcomm>& prcommList)
{

  int requestCount = 0;
  list<Prcomm>::iterator prcomm;
  for (prcomm = prcommList.begin(); prcomm!=prcommList.end(); prcomm++)
    if (prcomm->getNbrRank()!=mpi_rank) requestCount++;

  MPI_Request * sendRequestArray;
  MPI_Request * recvRequestArray;
  MPI_Status * statusArray;
  if (requestCount>0)
  {
    sendRequestArray = new MPI_Request[requestCount];
    recvRequestArray = new MPI_Request[requestCount];
    statusArray = new MPI_Status[requestCount];
  }

  // start...

  requestCount = 0;
  for (prcomm = prcommList.begin(); prcomm!=prcommList.end(); prcomm++)
  {

    // ensure size of packBufferDouble (we will be recving int pack in
    // a reverse exchange)...

    int npack = prcomm->packIndexVec.size();
    if (3*npack>prcomm->packBufferDoubleSize)
    {
      double * ptr;
      if ((ptr = (double*) realloc(prcomm->packBufferDouble, sizeof(double)*3*npack))==NULL)
      {
        cerr<<"Error: realloc failed in updateR2Reverse(double (*d)[3],..."<<endl;
        throw(-1);
      }
      prcomm->packBufferDouble = ptr;
      prcomm->packBufferDoubleSize = 3*npack;
    }

    // post the Irecv...

    if (prcomm->getNbrRank()!=mpi_rank)
    {
      MPI_Irecv(prcomm->packBufferDouble, 3*npack, MPI_DOUBLE, prcomm->getNbrRank(), UPDATE_R2_REVERSE_TAG, mpi_comm,
          &(recvRequestArray[requestCount]));
    }

    // ensure size of unpackBufferDouble...

    int nunpack = prcomm->unpackIndexVec.size();
    if (3*nunpack>prcomm->unpackBufferDoubleSize)
    {
      double * ptr;
      if ((ptr = (double*) realloc(prcomm->unpackBufferDouble, sizeof(double)*3*nunpack))==NULL)
      {
        cerr<<"Error: realloc failed in updateR2Reverse(double (*d)[3],..."<<endl;
        throw(-1);
      }
      prcomm->unpackBufferDouble = ptr;
      prcomm->unpackBufferDoubleSize = 3*nunpack;
    }

    // pack...
    switch (action)
    {
    case REPLACE_DATA:
    case ADD_DATA:
    case MIN_DATA:
    case MAX_DATA:
      for (int i = 0; i<nunpack; i++)
      {
        int index = prcomm->unpackIndexVec[i];
        for (int j = 0; j<3; j++)
          prcomm->unpackBufferDouble[3*i+j] = d[index][j];
      }
      break;
    case REPLACE_TRANSLATE_DATA:
    case ADD_TRANSLATE_DATA:
      for (list<Range>::iterator ri = prcomm->unpackRangeList.begin(); ri!=prcomm->unpackRangeList.end(); ri++)
      {
        for (int i = ri->getIndexFirst(); i<=ri->getIndexLast(); i++)
        {
          int index = prcomm->unpackIndexVec[i];
          for (int j = 0; j<3; j++)
            prcomm->unpackBufferDouble[3*i+j] = d[index][j];
        }
        // rangeTranslateVector is a virtual function that must be defined in 
        // any instantiated version of Gp...
        rangeTranslateVector(prcomm->unpackBufferDouble+3*ri->getIndexFirst(), // start
            ri->getIndexLast()-ri->getIndexFirst()+1, // size
            ri->getFlag(), ri->dxyz); // info about the transform
      }
      break;
    case REPLACE_ROTATE_DATA:
    case ADD_ROTATE_DATA:
      for (list<Range>::iterator ri = prcomm->unpackRangeList.begin(); ri!=prcomm->unpackRangeList.end(); ri++)
      {
        for (int i = ri->getIndexFirst(); i<=ri->getIndexLast(); i++)
        {
          int index = prcomm->unpackIndexVec[i];
          for (int j = 0; j<3; j++)
            prcomm->unpackBufferDouble[3*i+j] = d[index][j];
        }
        // rangeRotateVector is a virtual function that must be defined in 
        // any instantiated version of Gp...
        rangeRotateVector(prcomm->unpackBufferDouble+3*ri->getIndexFirst(), // start
            ri->getIndexLast()-ri->getIndexFirst()+1, // size
            ri->getFlag(), ri->dxyz); // info about the transform
      }
      break;
    default:
      cerr<<"Error: unsupported action: "<<action<<endl;
      throw(-1);
    }

    // and send...

    if (prcomm->getNbrRank()==mpi_rank)
    {
      assert(npack==nunpack);
      for (int i = 0; i<3*nunpack; i++)
        prcomm->packBufferDouble[i] = prcomm->unpackBufferDouble[i];
    }
    else
    {
      MPI_Issend(prcomm->unpackBufferDouble, 3*nunpack, MPI_DOUBLE, prcomm->getNbrRank(), UPDATE_R2_REVERSE_TAG,
          mpi_comm, &(sendRequestArray[requestCount]));
      // inc request count...
      requestCount++;
    }
  }

  // could do other work here - so eventually break up this routine into
  // 2 routines - start and finish...

  // now we wait for all messages to be received...

  if (requestCount>0) MPI_Waitall(requestCount, recvRequestArray, statusArray);

  // unpack...

  switch (action)
  {
  case REPLACE_DATA:
  case REPLACE_TRANSLATE_DATA:
  case REPLACE_ROTATE_DATA:
    for (prcomm = prcommList.begin(); prcomm!=prcommList.end(); prcomm++)
    {
      int npack = prcomm->packIndexVec.size();
      for (int i = 0; i<npack; i++)
      {
        int index = prcomm->packIndexVec[i];
        for (int j = 0; j<3; j++)
          d[index][j] = prcomm->packBufferDouble[3*i+j];
      }
    }
    break;
  case ADD_DATA:
  case ADD_TRANSLATE_DATA:
  case ADD_ROTATE_DATA:
    for (prcomm = prcommList.begin(); prcomm!=prcommList.end(); prcomm++)
    {
      int npack = prcomm->packIndexVec.size();
      for (int i = 0; i<npack; i++)
      {
        int index = prcomm->packIndexVec[i];
        for (int j = 0; j<3; j++)
          d[index][j] += prcomm->packBufferDouble[3*i+j];
      }
    }
    break;
  case MIN_DATA:
    for (prcomm = prcommList.begin(); prcomm!=prcommList.end(); prcomm++)
    {
      int npack = prcomm->packIndexVec.size();
      for (int i = 0; i<npack; i++)
      {
        int index = prcomm->packIndexVec[i];
        for (int j = 0; j<3; j++)
          d[index][j] = min(d[index][j], prcomm->packBufferDouble[3*i+j]);
      }
    }
    break;
  case MAX_DATA:
    for (prcomm = prcommList.begin(); prcomm!=prcommList.end(); prcomm++)
    {
      int npack = prcomm->packIndexVec.size();
      for (int i = 0; i<npack; i++)
      {
        int index = prcomm->packIndexVec[i];
        for (int j = 0; j<3; j++)
          d[index][j] = max(d[index][j], prcomm->packBufferDouble[3*i+j]);
      }
    }
    break;
  default:
    cerr<<"Error: unsupported action: "<<action<<endl;
    throw(-1);
  }

  // cleanup after sends...
  if (requestCount>0)
  {
    MPI_Waitall(requestCount, sendRequestArray, statusArray);
    delete[] sendRequestArray;
    delete[] recvRequestArray;
    delete[] statusArray;
  }

}

void Gp::updatePackBufferDoublePtrs(list<Prcomm>& prcommList)
{

  int requestCount = 0;
  list<Prcomm>::iterator prcomm;
  for (prcomm = prcommList.begin(); prcomm!=prcommList.end(); prcomm++)
    if (prcomm->getNbrRank()!=mpi_rank) requestCount++;

  MPI_Request * sendRequestArray;
  MPI_Request * recvRequestArray;
  MPI_Status * statusArray;
  if (requestCount>0)
  {
    sendRequestArray = new MPI_Request[requestCount];
    recvRequestArray = new MPI_Request[requestCount];
    statusArray = new MPI_Status[requestCount];
  }

  // start...

  requestCount = 0;
  for (prcomm = prcommList.begin(); prcomm!=prcommList.end(); prcomm++)
  {

    if (prcomm->getNbrRank()==mpi_rank)
    {
      prcomm->packBufferDoublePtr = prcomm->packBufferDouble;
    }
    else
    {
      MPI_Irecv(&prcomm->packBufferDoublePtr, 8, MPI_CHAR, prcomm->getNbrRank(), UPDATE_ADDRESS_TAG, mpi_comm,
          &(recvRequestArray[requestCount]));
      assert(prcomm->packBufferDoublePtr==NULL);
      cout<<"mpi_rank, nbr, address of packBufferDouble, test: "<<mpi_rank<<" "<<prcomm->getNbrRank()<<" "
          <<prcomm->packBufferDouble<<" "<<prcomm->packBufferDoublePtr<<endl;
      MPI_Issend(&prcomm->packBufferDouble, 8, MPI_CHAR, prcomm->getNbrRank(), UPDATE_ADDRESS_TAG, mpi_comm,
          &(sendRequestArray[requestCount]));
      // inc request count...
      requestCount++;
    }
  }

  // now we wait for all messages to be received...

  if (requestCount>0)
  {
    MPI_Waitall(requestCount, recvRequestArray, statusArray);
    delete[] sendRequestArray;
    delete[] recvRequestArray;
    delete[] statusArray;
  }

  for (prcomm = prcommList.begin(); prcomm!=prcommList.end(); prcomm++)
  {
    if (prcomm->getNbrRank()!=mpi_rank)
    {
      cout<<"mpi_rank, nbr, recved packBufferDoublePtr: "<<mpi_rank<<" "<<prcomm->getNbrRank()<<" "
          <<prcomm->packBufferDoublePtr<<endl;
    }
  }

  MPI_Pause("TAKE A LOOK");

}

void Gp::exchangePrcommBufferInt(list<Prcomm>& prcommList, const int n)
{

  assert(n>=1);

  int requestCount = 0;
  list<Prcomm>::iterator prcomm;
  for (prcomm = prcommList.begin(); prcomm!=prcommList.end(); prcomm++)
    if (prcomm->getNbrRank()!=mpi_rank) requestCount++;

  MPI_Request * sendRequestArray;
  MPI_Request * recvRequestArray;
  MPI_Status * statusArray;
  if (requestCount>0)
  {
    sendRequestArray = new MPI_Request[requestCount];
    recvRequestArray = new MPI_Request[requestCount];
    statusArray = new MPI_Status[requestCount];
  }

  // start...

  requestCount = 0;
  for (prcomm = prcommList.begin(); prcomm!=prcommList.end(); prcomm++)
  {

    // ensure size of unpackBufferInt...

    int nunpack = n*prcomm->unpackIndexVec.size();
    prcomm->ensureUnpackBufferIntSize(nunpack);

    // no packing neccessary - buffer is already packed...

    int npack = n*prcomm->packIndexVec.size();
    assert(npack<=prcomm->packBufferIntSize);

    // and send...

    if (prcomm->getNbrRank()==mpi_rank)
    {
      assert(npack==nunpack);
      int i;
      for (i = 0; i<npack; i++)
        prcomm->unpackBufferInt[i] = prcomm->packBufferInt[i];
    }
    else
    {
      MPI_Irecv(prcomm->unpackBufferInt, nunpack, MPI_INT, prcomm->getNbrRank(), EXCHANGE_INT_TAG, mpi_comm,
          &(recvRequestArray[requestCount]));
      MPI_Issend(prcomm->packBufferInt, npack, MPI_INT, prcomm->getNbrRank(), EXCHANGE_INT_TAG, mpi_comm,
          &(sendRequestArray[requestCount]));
      // inc request count...
      requestCount++;
    }
  }

  // now we wait for all messages to be received...

  if (requestCount>0)
  {
    MPI_Waitall(requestCount, recvRequestArray, statusArray);
    delete[] sendRequestArray;
    delete[] recvRequestArray;
    delete[] statusArray;
  }

}

void Gp::exchangePrcommBufferIntReverse(list<Prcomm>& prcommList, const int n)
{

  assert(n>=1);

  int requestCount = 0;
  list<Prcomm>::iterator prcomm;
  for (prcomm = prcommList.begin(); prcomm!=prcommList.end(); prcomm++)
    if (prcomm->getNbrRank()!=mpi_rank) requestCount++;

  MPI_Request * sendRequestArray;
  MPI_Request * recvRequestArray;
  MPI_Status * statusArray;
  if (requestCount>0)
  {
    sendRequestArray = new MPI_Request[requestCount];
    recvRequestArray = new MPI_Request[requestCount];
    statusArray = new MPI_Status[requestCount];
  }

  // start...

  requestCount = 0;
  for (prcomm = prcommList.begin(); prcomm!=prcommList.end(); prcomm++)
  {

    // ensure size of packBufferInt...

    int npack = n*prcomm->packIndexVec.size();
    prcomm->ensurePackBufferIntSize(npack);
    assert(npack>=0);
    // this is possible for asymmetric exchanges...
    //if ( npack == 0 )
    //  cout << "mpi_rank: " << mpi_rank << " prcomm->getNbrRank(): " << prcomm->getNbrRank() << " has npack == 0" << endl;

    // no packing neccessary - unpack buffer (recall this is a reverse
    // exchange) should already be packed...

    int nunpack = n*prcomm->unpackIndexVec.size();
    assert(nunpack>=0);
    assert(nunpack<=prcomm->unpackBufferIntSize);
    // this is possible for asymmetric exchanges...
    //if ( nunpack == 0 )
    //  cout << "mpi_rank: " << mpi_rank << " prcomm->getNbrRank(): " << prcomm->getNbrRank() << " has nunpack == 0" << endl;

    // and send...

    if (prcomm->getNbrRank()==mpi_rank)
    {
      assert(npack==nunpack);
      int i;
      for (i = 0; i<npack; i++)
        prcomm->packBufferInt[i] = prcomm->unpackBufferInt[i];
    }
    else
    {
      MPI_Irecv(prcomm->packBufferInt, npack, MPI_INT, prcomm->getNbrRank(), EXCHANGE_INT_REVERSE_TAG, mpi_comm,
          &(recvRequestArray[requestCount]));
      MPI_Issend(prcomm->unpackBufferInt, nunpack, MPI_INT, prcomm->getNbrRank(), EXCHANGE_INT_REVERSE_TAG, mpi_comm,
          &(sendRequestArray[requestCount]));
      // inc request count...
      requestCount++;
    }
  }

  // now we wait for all messages to be received...

  if (requestCount>0)
  {
    MPI_Waitall(requestCount, recvRequestArray, statusArray);
    delete[] sendRequestArray;
    delete[] recvRequestArray;
    delete[] statusArray;
  }

}

void Gp::updateI1(int * d, const int action, list<Prcomm>& prcommList)
{

  int requestCount = 0;
  list<Prcomm>::iterator prcomm;
  for (prcomm = prcommList.begin(); prcomm!=prcommList.end(); prcomm++)
    if (prcomm->getNbrRank()!=mpi_rank) requestCount++;

  MPI_Request * sendRequestArray;
  MPI_Request * recvRequestArray;
  MPI_Status * statusArray;
  if (requestCount>0)
  {
    sendRequestArray = new MPI_Request[requestCount];
    recvRequestArray = new MPI_Request[requestCount];
    statusArray = new MPI_Status[requestCount];
  }

  // start...

  requestCount = 0;
  for (prcomm = prcommList.begin(); prcomm!=prcommList.end(); prcomm++)
  {

    // ensure size of unpackBufferInt...

    int nunpack = prcomm->unpackIndexVec.size();
    prcomm->ensureUnpackBufferIntSize(nunpack);

    // post the Irecv...

    if (prcomm->getNbrRank()!=mpi_rank)
    {
      MPI_Irecv(prcomm->unpackBufferInt, nunpack, MPI_INT, prcomm->getNbrRank(), UPDATE_I1_TAG, mpi_comm,
          &(recvRequestArray[requestCount]));
      // XXX we should also store a pointer to the prcomm so when we use
      // MPI_Waitany in the future, we know exactly which prcomm is available
      // for unpacking.
    }

    // ensure size of packBufferInt...

    int npack = prcomm->packIndexVec.size();
    prcomm->ensurePackBufferIntSize(npack);

    // pack...

    int i;
    for (i = 0; i<npack; i++)
    {
      int index = prcomm->packIndexVec[i];
      prcomm->packBufferInt[i] = d[index];
    }

    // and send...

    if (prcomm->getNbrRank()==mpi_rank)
    {
      assert(npack==nunpack);
      for (i = 0; i<npack; i++)
        prcomm->unpackBufferInt[i] = prcomm->packBufferInt[i];
    }
    else
    {
      MPI_Issend(prcomm->packBufferInt, npack, MPI_INT, prcomm->getNbrRank(), UPDATE_I1_TAG, mpi_comm,
          &(sendRequestArray[requestCount]));
      // inc request count...
      requestCount++;
    }
  }

  // could do other work here - so eventually break up this routine into
  // 2 routines - start and finish...

  // now we wait for all messages to be received...

  if (requestCount>0) MPI_Waitall(requestCount, recvRequestArray, statusArray);

  // unpack...

  switch (action)
  {
  case REPLACE_DATA:
    for (prcomm = prcommList.begin(); prcomm!=prcommList.end(); prcomm++)
    {
      int nunpack = prcomm->unpackIndexVec.size();
      int i;
      for (i = 0; i<nunpack; i++)
      {
        int index = prcomm->unpackIndexVec[i];
        d[index] = prcomm->unpackBufferInt[i];
      }
    }
    break;
  case ADD_DATA:
    for (prcomm = prcommList.begin(); prcomm!=prcommList.end(); prcomm++)
    {
      int nunpack = prcomm->unpackIndexVec.size();
      int i;
      for (i = 0; i<nunpack; i++)
      {
        int index = prcomm->unpackIndexVec[i];
        d[index] += prcomm->unpackBufferInt[i];
      }
    }
    break;
  case MIN_DATA:
    for (prcomm = prcommList.begin(); prcomm!=prcommList.end(); prcomm++)
    {
      int nunpack = prcomm->unpackIndexVec.size();
      int i;
      for (i = 0; i<nunpack; i++)
      {
        int index = prcomm->unpackIndexVec[i];
        d[index] = min(d[index], prcomm->unpackBufferInt[i]);
      }
    }
    break;
  case MAX_DATA:
    for (prcomm = prcommList.begin(); prcomm!=prcommList.end(); prcomm++)
    {
      int nunpack = prcomm->unpackIndexVec.size();
      int i;
      for (i = 0; i<nunpack; i++)
      {
        int index = prcomm->unpackIndexVec[i];
        d[index] = max(d[index], prcomm->unpackBufferInt[i]);
      }
    }
    break;
  case MIN_NO_PERIODIC_DATA:
    for (prcomm = prcommList.begin(); prcomm!=prcommList.end(); prcomm++)
    {
      list<Range>::iterator range;
      for (range = prcomm->unpackRangeList.begin(); range!=prcomm->unpackRangeList.end(); range++)
      {
        // non-periodic ranges only...
        if (range->isPeriodic()==0)
        {
          int i;
          for (i = range->getIndexFirst(); i<=range->getIndexLast(); i++)
          {
            int index = prcomm->unpackIndexVec[i];
            d[index] = min(d[index], prcomm->unpackBufferInt[i]);
          }
        }
      }
    }
    break;
  case MAX_NO_PERIODIC_DATA:
    for (prcomm = prcommList.begin(); prcomm!=prcommList.end(); prcomm++)
    {
      list<Range>::iterator range;
      for (range = prcomm->unpackRangeList.begin(); range!=prcomm->unpackRangeList.end(); range++)
      {
        // non-periodic ranges only...
        if (range->isPeriodic()==0)
        {
          int i;
          for (i = range->getIndexFirst(); i<=range->getIndexLast(); i++)
          {
            int index = prcomm->unpackIndexVec[i];
            d[index] = max(d[index], prcomm->unpackBufferInt[i]);
          }
        }
      }
    }
    break;
  case ADD_NO_PERIODIC_DATA:
    for (prcomm = prcommList.begin(); prcomm!=prcommList.end(); prcomm++)
    {
      list<Range>::iterator range;
      for (range = prcomm->unpackRangeList.begin(); range!=prcomm->unpackRangeList.end(); range++)
      {
        // non-periodic ranges only...
        if (range->isPeriodic()==0)
        {
          int i;
          for (i = range->getIndexFirst(); i<=range->getIndexLast(); i++)
          {
            int index = prcomm->unpackIndexVec[i];
            d[index] += prcomm->unpackBufferInt[i];
          }
        }
      }
    }
    break;
  case BITWISE_OR_DATA:
    for (prcomm = prcommList.begin(); prcomm!=prcommList.end(); prcomm++)
    {
      int nunpack = prcomm->unpackIndexVec.size();
      int i;
      for (i = 0; i<nunpack; i++)
      {
        int index = prcomm->unpackIndexVec[i];
        d[index] |= prcomm->unpackBufferInt[i];
      }
    }
    break;
  case NO_CHANGE_DATA:
    // fills the buffers, but makes no changes...
    break;
  default:
    cerr<<"Error: unsupported action: "<<action<<endl;
    throw(-1);
  }

  // cleanup after sends...
  if (requestCount>0)
  {
    MPI_Waitall(requestCount, sendRequestArray, statusArray);
    delete[] sendRequestArray;
    delete[] recvRequestArray;
    delete[] statusArray;
  }

}

void Gp::updateI1Reverse(int * d, const int action, list<Prcomm>& prcommList)
{

  int requestCount = 0;
  list<Prcomm>::iterator prcomm;
  for (prcomm = prcommList.begin(); prcomm!=prcommList.end(); prcomm++)
    if (prcomm->getNbrRank()!=mpi_rank) requestCount++;

  MPI_Request * sendRequestArray;
  MPI_Request * recvRequestArray;
  MPI_Status * statusArray;
  if (requestCount>0)
  {
    sendRequestArray = new MPI_Request[requestCount];
    recvRequestArray = new MPI_Request[requestCount];
    statusArray = new MPI_Status[requestCount];
  }

  // start...

  requestCount = 0;
  for (prcomm = prcommList.begin(); prcomm!=prcommList.end(); prcomm++)
  {

    // ensure size of packBufferInt (used for unpacking in a reverse exchange)...

    int nunpack = prcomm->packIndexVec.size();
    if (nunpack>prcomm->packBufferIntSize)
    {
      int * ptr;
      if ((ptr = (int*) realloc(prcomm->packBufferInt, sizeof(int)*nunpack))==NULL)
      {
        cerr<<"Error: realloc failed in updateI1(int * d,..."<<endl;
        throw(-1);
      }
      prcomm->packBufferInt = ptr;
      prcomm->packBufferIntSize = nunpack;
    }

    // post the Irecv...

    if (prcomm->getNbrRank()!=mpi_rank)
    {
      MPI_Irecv(prcomm->packBufferInt, nunpack, MPI_INT, prcomm->getNbrRank(), UPDATE_I1_REVERSE_TAG, mpi_comm,
          &(recvRequestArray[requestCount]));
      // XXX we should also store a pointer to the prcomm so when we use
      // MPI_Waitany in the future, we know exactly which prcomm is available
      // for unpacking.
    }

    // ensure size of unpackBufferInt...

    int npack = prcomm->unpackIndexVec.size();
    if (npack>prcomm->unpackBufferIntSize)
    {
      int * ptr;
      if ((ptr = (int*) realloc(prcomm->unpackBufferInt, sizeof(int)*npack))==NULL)
      {
        cerr<<"Error: realloc failed in updateI1"<<endl;
        throw(-1);
      }
      prcomm->unpackBufferInt = ptr;
      prcomm->unpackBufferIntSize = npack;
    }

    // pack...

    int i;
    for (i = 0; i<npack; i++)
    {
      int index = prcomm->unpackIndexVec[i];
      prcomm->unpackBufferInt[i] = d[index];
    }

    // and send...

    if (prcomm->getNbrRank()==mpi_rank)
    {
      assert(npack==nunpack);
      for (i = 0; i<npack; i++)
        prcomm->packBufferInt[i] = prcomm->unpackBufferInt[i];
    }
    else
    {
      MPI_Issend(prcomm->unpackBufferInt, npack, MPI_INT, prcomm->getNbrRank(), UPDATE_I1_REVERSE_TAG, mpi_comm,
          &(sendRequestArray[requestCount]));
      // inc request count...
      requestCount++;
    }
  }

  // ---------------------------------------------------------------------
  // could do other work here - so eventually break up this routine into
  // 2 routines - start and finish...
  // ---------------------------------------------------------------------

  // now we wait for all messages to be received...

  if (requestCount>0) MPI_Waitall(requestCount, recvRequestArray, statusArray);

  // unpack...

  switch (action)
  {
  case REPLACE_DATA:
    for (prcomm = prcommList.begin(); prcomm!=prcommList.end(); prcomm++)
    {
      int nunpack = prcomm->packIndexVec.size();
      int i;
      for (i = 0; i<nunpack; i++)
      {
        int index = prcomm->packIndexVec[i];
        d[index] = prcomm->packBufferInt[i];
      }
    }
    break;
  case ADD_DATA:
    for (prcomm = prcommList.begin(); prcomm!=prcommList.end(); prcomm++)
    {
      int nunpack = prcomm->packIndexVec.size();
      int i;
      for (i = 0; i<nunpack; i++)
      {
        int index = prcomm->packIndexVec[i];
        d[index] += prcomm->packBufferInt[i];
      }
    }
    break;
  case MIN_DATA:
    for (prcomm = prcommList.begin(); prcomm!=prcommList.end(); prcomm++)
    {
      int nunpack = prcomm->packIndexVec.size();
      int i;
      for (i = 0; i<nunpack; i++)
      {
        int index = prcomm->packIndexVec[i];
        d[index] = min(d[index], prcomm->packBufferInt[i]);
      }
    }
    break;
  case MAX_DATA:
    for (prcomm = prcommList.begin(); prcomm!=prcommList.end(); prcomm++)
    {
      int nunpack = prcomm->packIndexVec.size();
      int i;
      for (i = 0; i<nunpack; i++)
      {
        int index = prcomm->packIndexVec[i];
        d[index] = max(d[index], prcomm->packBufferInt[i]);
      }
    }
    break;
  case BITWISE_OR_DATA:
    for (prcomm = prcommList.begin(); prcomm!=prcommList.end(); prcomm++)
    {
      int nunpack = prcomm->packIndexVec.size();
      int i;
      for (i = 0; i<nunpack; i++)
      {
        int index = prcomm->packIndexVec[i];
        d[index] |= prcomm->packBufferInt[i];
      }
    }
    break;
  case NO_CHANGE_DATA:
    // fills the buffers, but makes no changes...
    break;
  default:
    cerr<<"Error: unsupported action: "<<action<<endl;
    throw(-1);
  }

  // cleanup after sends...
  if (requestCount>0)
  {
    MPI_Waitall(requestCount, sendRequestArray, statusArray);
    delete[] sendRequestArray;
    delete[] recvRequestArray;
    delete[] statusArray;
  }

}

void Gp::updateI2(int(*d)[3], const int action, list<Prcomm>& prcommList)
{

  // count the requests, not including self-requests, which do not 
  // invoke the MPI communication. Requests have to be put into an array
  // ..

  int requestCount = 0;
  list<Prcomm>::iterator prcomm;
  for (prcomm = prcommList.begin(); prcomm!=prcommList.end(); prcomm++)
    if (prcomm->getNbrRank()!=mpi_rank) requestCount++;

  MPI_Request * sendRequestArray;
  MPI_Request * recvRequestArray;
  MPI_Status * statusArray;
  if (requestCount>0)
  {
    sendRequestArray = new MPI_Request[requestCount];
    recvRequestArray = new MPI_Request[requestCount];
    statusArray = new MPI_Status[requestCount];
  }

  // start...

  requestCount = 0;
  for (prcomm = prcommList.begin(); prcomm!=prcommList.end(); prcomm++)
  {

    // ensure size of unpackBufferInt...

    int nunpack = prcomm->unpackIndexVec.size();
    if (3*nunpack>prcomm->unpackBufferIntSize)
    {
      int * ptr;
      if ((ptr = (int*) realloc(prcomm->unpackBufferInt, sizeof(int)*3*nunpack))==NULL)
      {
        cerr<<"Error: realloc failed in updateI2(int (*d)[3],..."<<endl;
        throw(-1);
      }
      prcomm->unpackBufferInt = ptr;
      prcomm->unpackBufferIntSize = 3*nunpack;
    }

    // post the Irecv...

    if (prcomm->getNbrRank()!=mpi_rank)
    {
      MPI_Irecv(prcomm->unpackBufferInt, 3*nunpack, MPI_INT, prcomm->getNbrRank(), UPDATE_I2_TAG, mpi_comm,
          &(recvRequestArray[requestCount]));
      // XXX we should also store a pointer to the prcomm so when we use
      // MPI_Waitany in the future, we know exactly which prcomm is available
      // for unpacking.
    }

    // ensure size of packBufferInt...

    int npack = prcomm->packIndexVec.size();
    if (3*npack>prcomm->packBufferIntSize)
    {
      int * ptr;
      if ((ptr = (int*) realloc(prcomm->packBufferInt, sizeof(int)*3*npack))==NULL)
      {
        cerr<<"Error: realloc failed in updateI2"<<endl;
        throw(-1);
      }
      prcomm->packBufferInt = ptr;
      prcomm->packBufferIntSize = 3*npack;
    }

    // pack...

    int i;
    for (i = 0; i<npack; i++)
    {
      int index = prcomm->packIndexVec[i];
      int j;
      for (j = 0; j<3; j++)
        prcomm->packBufferInt[3*i+j] = d[index][j];
    }

    // and send...

    if (prcomm->getNbrRank()==mpi_rank)
    {
      assert(npack==nunpack);
      for (i = 0; i<3*npack; i++)
        prcomm->unpackBufferInt[i] = prcomm->packBufferInt[i];
    }
    else
    {
      MPI_Issend(prcomm->packBufferInt, 3*npack, MPI_INT, prcomm->getNbrRank(), UPDATE_I2_TAG, mpi_comm,
          &(sendRequestArray[requestCount]));
      // inc request count...
      requestCount++;
    }
  }

  // could do other work here - so eventually break up this routine into
  // 2 routines - start and finish...

  // now we wait for all messages to be received...

  if (requestCount>0) MPI_Waitall(requestCount, recvRequestArray, statusArray);

  // unpack...

  switch (action)
  {
  case REPLACE_DATA:
    for (prcomm = prcommList.begin(); prcomm!=prcommList.end(); prcomm++)
    {
      int nunpack = prcomm->unpackIndexVec.size();
      int i;
      for (i = 0; i<nunpack; i++)
      {
        int index = prcomm->unpackIndexVec[i];
        int j;
        for (j = 0; j<3; j++)
          d[index][j] = prcomm->unpackBufferInt[3*i+j];
      }
    }
    break;
  case ADD_DATA:
    for (prcomm = prcommList.begin(); prcomm!=prcommList.end(); prcomm++)
    {
      int nunpack = prcomm->unpackIndexVec.size();
      int i;
      for (i = 0; i<nunpack; i++)
      {
        int index = prcomm->unpackIndexVec[i];
        int j;
        for (j = 0; j<3; j++)
          d[index][j] += prcomm->unpackBufferInt[3*i+j];
      }
    }
    break;
  case MIN_DATA:
    for (prcomm = prcommList.begin(); prcomm!=prcommList.end(); prcomm++)
    {
      int nunpack = prcomm->unpackIndexVec.size();
      int i;
      for (i = 0; i<nunpack; i++)
      {
        int index = prcomm->unpackIndexVec[i];
        int j;
        for (j = 0; j<3; j++)
          d[index][j] = min(d[index][j], prcomm->unpackBufferInt[3*i+j]);
      }
    }
    break;
  case MAX_DATA:
    for (prcomm = prcommList.begin(); prcomm!=prcommList.end(); prcomm++)
    {
      int nunpack = prcomm->unpackIndexVec.size();
      int i;
      for (i = 0; i<nunpack; i++)
      {
        int index = prcomm->unpackIndexVec[i];
        int j;
        for (j = 0; j<3; j++)
          d[index][j] = max(d[index][j], prcomm->unpackBufferInt[3*i+j]);
      }
    }
    break;
  default:
    cerr<<"Error: unsupported action: "<<action<<endl;
    throw(-1);
  }

  // cleanup after sends...
  if (requestCount>0)
  {
    MPI_Waitall(requestCount, sendRequestArray, statusArray);
    delete[] sendRequestArray;
    delete[] recvRequestArray;
    delete[] statusArray;
  }

}

void Gp::updateI2Reverse(int(*d)[3], const int action, list<Prcomm>& prcommList)
{

  // count the requests, not including self-requests, which do not 
  // invoke the MPI communication. Requests have to be put into an array
  // ..

  int requestCount = 0;
  list<Prcomm>::iterator prcomm;
  for (prcomm = prcommList.begin(); prcomm!=prcommList.end(); prcomm++)
    if (prcomm->getNbrRank()!=mpi_rank) requestCount++;

  MPI_Request * sendRequestArray;
  MPI_Request * recvRequestArray;
  MPI_Status * statusArray;
  if (requestCount>0)
  {
    sendRequestArray = new MPI_Request[requestCount];
    recvRequestArray = new MPI_Request[requestCount];
    statusArray = new MPI_Status[requestCount];
  }

  // start...

  requestCount = 0;
  for (prcomm = prcommList.begin(); prcomm!=prcommList.end(); prcomm++)
  {

    // ensure size of packBufferInt (used for unpacking in these reverse routines)...

    int nunpack = prcomm->packIndexVec.size();
    if (3*nunpack>prcomm->packBufferIntSize)
    {
      int * ptr;
      if ((ptr = (int*) realloc(prcomm->packBufferInt, sizeof(int)*3*nunpack))==NULL)
      {
        cerr<<"Error: realloc failed in updateI2Reverse(int (*d)[3],..."<<endl;
        throw(-1);
      }
      prcomm->packBufferInt = ptr;
      prcomm->packBufferIntSize = 3*nunpack;
    }

    // post the Irecv...

    if (prcomm->getNbrRank()!=mpi_rank)
    {
      MPI_Irecv(prcomm->packBufferInt, 3*nunpack, MPI_INT, prcomm->getNbrRank(), UPDATE_I2_REVERSE_TAG, mpi_comm,
          &(recvRequestArray[requestCount]));
      // XXX we should also store a pointer to the prcomm so when we use
      // MPI_Waitany in the future, we know exactly which prcomm is available
      // for unpacking.
    }

    // ensure size of unpackBufferInt...

    int npack = prcomm->unpackIndexVec.size();
    if (3*npack>prcomm->unpackBufferIntSize)
    {
      int * ptr;
      if ((ptr = (int*) realloc(prcomm->unpackBufferInt, sizeof(int)*3*npack))==NULL)
      {
        cerr<<"Error: realloc failed in updateI2Reverse"<<endl;
        throw(-1);
      }
      prcomm->unpackBufferInt = ptr;
      prcomm->unpackBufferIntSize = 3*npack;
    }

    // pack...

    int i;
    for (i = 0; i<npack; i++)
    {
      int index = prcomm->unpackIndexVec[i];
      int j;
      for (j = 0; j<3; j++)
        prcomm->unpackBufferInt[3*i+j] = d[index][j];
    }

    // and send...

    if (prcomm->getNbrRank()==mpi_rank)
    {
      assert(npack==nunpack);
      for (i = 0; i<3*npack; i++)
        prcomm->packBufferInt[i] = prcomm->unpackBufferInt[i];
    }
    else
    {
      MPI_Issend(prcomm->unpackBufferInt, 3*npack, MPI_INT, prcomm->getNbrRank(), UPDATE_I2_REVERSE_TAG, mpi_comm,
          &(sendRequestArray[requestCount]));
      // inc request count...
      requestCount++;
    }
  }

  // could do other work here - so eventually break up this routine into
  // 2 routines - start and finish...

  // now we wait for all messages to be received...

  if (requestCount>0) MPI_Waitall(requestCount, recvRequestArray, statusArray);

  // unpack...

  switch (action)
  {
  case REPLACE_DATA:
    for (prcomm = prcommList.begin(); prcomm!=prcommList.end(); prcomm++)
    {
      int nunpack = prcomm->packIndexVec.size();
      int i;
      for (i = 0; i<nunpack; i++)
      {
        int index = prcomm->packIndexVec[i];
        int j;
        for (j = 0; j<3; j++)
          d[index][j] = prcomm->packBufferInt[3*i+j];
      }
    }
    break;
  case ADD_DATA:
    for (prcomm = prcommList.begin(); prcomm!=prcommList.end(); prcomm++)
    {
      int nunpack = prcomm->packIndexVec.size();
      int i;
      for (i = 0; i<nunpack; i++)
      {
        int index = prcomm->packIndexVec[i];
        int j;
        for (j = 0; j<3; j++)
          d[index][j] += prcomm->packBufferInt[3*i+j];
      }
    }
    break;
  case MIN_DATA:
    for (prcomm = prcommList.begin(); prcomm!=prcommList.end(); prcomm++)
    {
      int nunpack = prcomm->packIndexVec.size();
      int i;
      for (i = 0; i<nunpack; i++)
      {
        int index = prcomm->packIndexVec[i];
        int j;
        for (j = 0; j<3; j++)
          d[index][j] = min(d[index][j], prcomm->packBufferInt[3*i+j]);
      }
    }
    break;
  case MAX_DATA:
    for (prcomm = prcommList.begin(); prcomm!=prcommList.end(); prcomm++)
    {
      int nunpack = prcomm->packIndexVec.size();
      int i;
      for (i = 0; i<nunpack; i++)
      {
        int index = prcomm->packIndexVec[i];
        int j;
        for (j = 0; j<3; j++)
          d[index][j] = max(d[index][j], prcomm->packBufferInt[3*i+j]);
      }
    }
    break;
  default:
    cerr<<"Error: unsupported action: "<<action<<endl;
    throw(-1);
  }

  // cleanup after sends...
  if (requestCount>0)
  {
    MPI_Waitall(requestCount, sendRequestArray, statusArray);
    delete[] sendRequestArray;
    delete[] recvRequestArray;
    delete[] statusArray;
  }

}

void Gp::exchangePrcommNpackV(list<Prcomm>& prcommList)
{

  // exchange the pack_v to the corresponding unpack_v on paired
  // communicators of prcommmList...

  int requestCount = 0;
  list<Prcomm>::iterator prcomm;
  for (prcomm = prcommList.begin(); prcomm!=prcommList.end(); prcomm++)
    if (prcomm->getNbrRank()!=mpi_rank) requestCount++;

  MPI_Request * sendRequestArray;
  MPI_Request * recvRequestArray;
  MPI_Status * statusArray;
  if (requestCount>0)
  {
    sendRequestArray = new MPI_Request[requestCount];
    recvRequestArray = new MPI_Request[requestCount];
    statusArray = new MPI_Status[requestCount];
  }

  // start...

  requestCount = 0;
  for (prcomm = prcommList.begin(); prcomm!=prcommList.end(); prcomm++)
  {

    // post the Irecv...
    if (prcomm->getNbrRank()!=mpi_rank)
    {
      MPI_Irecv(&(prcomm->nunpack_v), 1, MPI_INT, prcomm->getNbrRank(), EXCHANGE_PRCOMM_NPACKV_TAG, mpi_comm,
          &(recvRequestArray[requestCount]));
    }

    // and send...      
    if (prcomm->getNbrRank()==mpi_rank)
    {
      prcomm->nunpack_v = prcomm->npack_v;
    }
    else
    {
      MPI_Issend(&(prcomm->npack_v), 1, MPI_INT, prcomm->getNbrRank(), EXCHANGE_PRCOMM_NPACKV_TAG, mpi_comm,
          &(sendRequestArray[requestCount]));
      // inc request count...
      requestCount++;
    }
  }

  // now we wait for all messages to be received (and sent)...
  if (requestCount>0)
  {
    MPI_Waitall(requestCount, recvRequestArray, statusArray);
    MPI_Waitall(requestCount, sendRequestArray, statusArray);
    delete[] sendRequestArray;
    delete[] recvRequestArray;
    delete[] statusArray;
  }
}

void Gp::exchangePrcommNunpackV(list<Prcomm>& prcommList)
{

  // exchange the nunpack_v to the corresponding npack_v on paired
  // communicators of prcommmList...

  int requestCount = 0;
  list<Prcomm>::iterator prcomm;
  for (prcomm = prcommList.begin(); prcomm!=prcommList.end(); prcomm++)
    if (prcomm->getNbrRank()!=mpi_rank) requestCount++;

  MPI_Request * sendRequestArray;
  MPI_Request * recvRequestArray;
  MPI_Status * statusArray;
  if (requestCount>0)
  {
    sendRequestArray = new MPI_Request[requestCount];
    recvRequestArray = new MPI_Request[requestCount];
    statusArray = new MPI_Status[requestCount];
  }

  // start...

  requestCount = 0;
  for (prcomm = prcommList.begin(); prcomm!=prcommList.end(); prcomm++)
  {

    // post the Irecv...
    if (prcomm->getNbrRank()!=mpi_rank)
    {
      MPI_Irecv(&(prcomm->npack_v), 1, MPI_INT, prcomm->getNbrRank(), EXCHANGE_PRCOMM_NUNPACKV_TAG, mpi_comm,
          &(recvRequestArray[requestCount]));
    }

    // and send...      
    if (prcomm->getNbrRank()==mpi_rank)
    {
      prcomm->npack_v = prcomm->nunpack_v;
    }
    else
    {
      MPI_Issend(&(prcomm->nunpack_v), 1, MPI_INT, prcomm->getNbrRank(), EXCHANGE_PRCOMM_NUNPACKV_TAG, mpi_comm,
          &(sendRequestArray[requestCount]));
      // inc request count...
      requestCount++;
    }
  }

  // now we wait for all messages to be received (and sent)...
  if (requestCount>0)
  {
    MPI_Waitall(requestCount, recvRequestArray, statusArray);
    MPI_Waitall(requestCount, sendRequestArray, statusArray);
    delete[] sendRequestArray;
    delete[] recvRequestArray;
    delete[] statusArray;
  }
}

void Gp::exchangePrcommBufferIntV(list<Prcomm>& prcommList, const int n)
{

  // the passed prcommList has had its packBufferInt's packed with npack_v*n ints, 
  // so update with nbrs...

  int requestCount = 0;
  list<Prcomm>::iterator prcomm;
  for (prcomm = prcommList.begin(); prcomm!=prcommList.end(); prcomm++)
    if (prcomm->getNbrRank()!=mpi_rank) requestCount++;

  MPI_Request * sendRequestArray;
  MPI_Request * recvRequestArray;
  MPI_Status * statusArray;
  if (requestCount>0)
  {
    sendRequestArray = new MPI_Request[requestCount];
    recvRequestArray = new MPI_Request[requestCount];
    statusArray = new MPI_Status[requestCount];
  }

  // start...

  requestCount = 0;
  for (prcomm = prcommList.begin(); prcomm!=prcommList.end(); prcomm++)
  {

    // post the Irecv...
    if (prcomm->getNbrRank()!=mpi_rank)
    {
      MPI_Irecv(&(prcomm->nunpack_v), 1, MPI_INT, prcomm->getNbrRank(), EXCHANGE_PRCOMM_INTV_TAG1, mpi_comm,
          &(recvRequestArray[requestCount]));
    }

    // and send...      
    if (prcomm->getNbrRank()==mpi_rank)
    {
      prcomm->nunpack_v = prcomm->npack_v;
    }
    else
    {
      MPI_Issend(&(prcomm->npack_v), 1, MPI_INT, prcomm->getNbrRank(), EXCHANGE_PRCOMM_INTV_TAG1, mpi_comm,
          &(sendRequestArray[requestCount]));
      // inc request count...
      requestCount++;
    }
  }

  // now we wait for all messages to be received (and sent)...
  if (requestCount>0)
  {
    MPI_Waitall(requestCount, recvRequestArray, statusArray);
    MPI_Waitall(requestCount, sendRequestArray, statusArray);
  }

  // now each recv knows how much data to expect, so make a new 
  // send and recv requestCount's so we don't send any message associated
  // with zero-sized data...

  int sendRequestCount = 0;
  int recvRequestCount = 0;
  for (prcomm = prcommList.begin(); prcomm!=prcommList.end(); prcomm++)
  {

    // post the Irecv...
    if (prcomm->nunpack_v>0)
    {
      prcomm->ensureUnpackBufferIntSize(prcomm->nunpack_v*n);
      if (prcomm->getNbrRank()!=mpi_rank)
      {
        MPI_Irecv(prcomm->unpackBufferInt, prcomm->nunpack_v*n, MPI_INT, prcomm->getNbrRank(),
            EXCHANGE_PRCOMM_INTV_TAG2, mpi_comm, &(recvRequestArray[recvRequestCount]));
        recvRequestCount++;
      }
    }

    // and send...      
    if (prcomm->npack_v>0)
    {
      if (prcomm->getNbrRank()==mpi_rank)
      {
        int i;
        for (i = 0; i<prcomm->npack_v*n; i++)
          prcomm->unpackBufferInt[i] = prcomm->packBufferInt[i];
      }
      else
      {
        MPI_Issend(prcomm->packBufferInt, prcomm->npack_v*n, MPI_INT, prcomm->getNbrRank(), EXCHANGE_PRCOMM_INTV_TAG2,
            mpi_comm, &(sendRequestArray[sendRequestCount]));
        // inc request count...
        sendRequestCount++;
      }
    }
  }

  if (requestCount>0)
  {
    if (recvRequestCount>0) MPI_Waitall(recvRequestCount, recvRequestArray, statusArray);
    if (sendRequestCount>0) MPI_Waitall(sendRequestCount, sendRequestArray, statusArray);
    delete[] sendRequestArray;
    delete[] recvRequestArray;
    delete[] statusArray;
  }

}

void Gp::exchangePrcommBufferDoubleV(list<Prcomm>& prcommList, const int n)
{

  // the passed prcommList has had its packBufferDouble's packed with npack_v*n doubles, 
  // so update with nbrs...

  int requestCount = 0;
  list<Prcomm>::iterator prcomm;
  for (prcomm = prcommList.begin(); prcomm!=prcommList.end(); prcomm++)
    if (prcomm->getNbrRank()!=mpi_rank) requestCount++;

  MPI_Request * sendRequestArray;
  MPI_Request * recvRequestArray;
  MPI_Status * statusArray;
  if (requestCount>0)
  {
    sendRequestArray = new MPI_Request[requestCount];
    recvRequestArray = new MPI_Request[requestCount];
    statusArray = new MPI_Status[requestCount];
  }

  // start...

  requestCount = 0;
  for (prcomm = prcommList.begin(); prcomm!=prcommList.end(); prcomm++)
  {

    // post the Irecv...
    if (prcomm->getNbrRank()!=mpi_rank)
    {
      MPI_Irecv(&(prcomm->nunpack_v), 1, MPI_INT, prcomm->getNbrRank(), EXCHANGE_PRCOMM_DOUBLEV_TAG1, mpi_comm,
          &(recvRequestArray[requestCount]));
    }

    // and send...      
    if (prcomm->getNbrRank()==mpi_rank)
    {
      prcomm->nunpack_v = prcomm->npack_v;
    }
    else
    {
      MPI_Issend(&(prcomm->npack_v), 1, MPI_INT, prcomm->getNbrRank(), EXCHANGE_PRCOMM_DOUBLEV_TAG1, mpi_comm,
          &(sendRequestArray[requestCount]));
      // inc request count...
      requestCount++;
    }
  }

  // now we wait for all messages to be received (and sent)...
  if (requestCount>0)
  {
    MPI_Waitall(requestCount, recvRequestArray, statusArray);
    MPI_Waitall(requestCount, sendRequestArray, statusArray);
  }

  // now each recv knows how much data to expect, so make a new 
  // send and recv requestCount's so we don't send any message associated
  // with zero-sized data...

  int sendRequestCount = 0;
  int recvRequestCount = 0;
  for (prcomm = prcommList.begin(); prcomm!=prcommList.end(); prcomm++)
  {

    // post the Irecv...
    if (prcomm->nunpack_v>0)
    {
      prcomm->ensureUnpackBufferDoubleSize(prcomm->nunpack_v*n);
      if (prcomm->getNbrRank()!=mpi_rank)
      {
        MPI_Irecv(prcomm->unpackBufferDouble, prcomm->nunpack_v*n, MPI_DOUBLE, prcomm->getNbrRank(),
            EXCHANGE_PRCOMM_DOUBLEV_TAG2, mpi_comm, &(recvRequestArray[recvRequestCount]));
        recvRequestCount++;
      }
    }

    // and send...      
    if (prcomm->npack_v>0)
    {
      if (prcomm->getNbrRank()==mpi_rank)
      {
        for (int i = 0; i<prcomm->npack_v*n; i++)
          prcomm->unpackBufferDouble[i] = prcomm->packBufferDouble[i];
      }
      else
      {
        MPI_Issend(prcomm->packBufferDouble, prcomm->npack_v*n, MPI_DOUBLE, prcomm->getNbrRank(),
            EXCHANGE_PRCOMM_DOUBLEV_TAG2, mpi_comm, &(sendRequestArray[sendRequestCount]));
        // inc request count...
        sendRequestCount++;
      }
    }
  }

  if (requestCount>0)
  {
    if (recvRequestCount>0) MPI_Waitall(recvRequestCount, recvRequestArray, statusArray);
    if (sendRequestCount>0) MPI_Waitall(sendRequestCount, sendRequestArray, statusArray);
    delete[] sendRequestArray;
    delete[] recvRequestArray;
    delete[] statusArray;
  }

}

void Gp::dumpFacePrcommData()
{

  int rank;
  for (rank = 0; rank<mpi_size; rank++)
  {

    if (mpi_rank==rank)
    {

      cout<<mpi_rank<<":dumpFacePrcommData()..."<<endl;
      cout<<mpi_rank<<":facePrcommList.size(): "<<facePrcommList.size()<<endl;
      cout<<mpi_rank<<":face prcomms:"<<endl;
      list<Prcomm>::iterator prcomm;
      for (prcomm = facePrcommList.begin(); prcomm!=facePrcommList.end(); prcomm++)
      {
        cout<<mpi_rank<<":  =================================== "<<endl;
        cout<<mpi_rank<<":  prcomm->getNbrRank()           : "<<prcomm->getNbrRank()<<endl;
        cout<<mpi_rank<<":  prcomm->packIndexVec.size()    : "<<prcomm->packIndexVec.size()<<endl;
        cout<<mpi_rank<<":  prcomm->packIndexVec           : ";
        int i;
        for (i = 0; i<min(12, (int) prcomm->packIndexVec.size()); i++)
          cout<<prcomm->packIndexVec[i]<<" ";
        if (prcomm->packIndexVec.size()>12) cout<<mpi_rank<<":..."<<endl;
        else cout<<endl;
        cout<<mpi_rank<<":  prcomm->packRangeList.size()   : "<<prcomm->packRangeList.size()<<endl;
        cout<<mpi_rank<<":  pack ranges:"<<endl;
        int i_f = 0;
        list<Range>::iterator range;
        for (range = prcomm->packRangeList.begin(); range!=prcomm->packRangeList.end(); range++)
        {
          cout<<mpi_rank<<":    -------------------------------- "<<endl;
          cout<<mpi_rank<<":    range->getBits()             : ";
          int ii;
          for (ii = 0; ii<32; ii++)
          {
            if (range->getBits()&(1<<ii))
            {
              cout<<"1";
            }
            else
            {
              cout<<"0";
            }
          }
          cout<<endl;
          cout<<mpi_rank<<":    range->getIndexFirst/Last    : "<<range->getIndexFirst()<<" "<<range->getIndexLast()
              <<endl;
          assert(range->getIndexFirst()==i_f);
          i_f = range->getIndexLast()+1;
        }
        assert(prcomm->packIndexVec.size()==i_f);
        cout<<mpi_rank<<":    -------------------------------- "<<endl;
        cout<<mpi_rank<<":  end of pack ranges"<<endl;
        cout<<mpi_rank<<":  prcomm->unpackRangeList.size() : "<<prcomm->unpackRangeList.size()<<endl;
        cout<<mpi_rank<<":  prcomm->unpackIndexVec.size()  : "<<prcomm->unpackIndexVec.size()<<endl;
        cout<<mpi_rank<<":  prcomm->unpackIndexVec         : ";
        for (i = 0; i<min(12, (int) prcomm->unpackIndexVec.size()); i++)
          cout<<prcomm->unpackIndexVec[i]<<" ";
        if (prcomm->unpackIndexVec.size()>12) cout<<mpi_rank<<":..."<<endl;
        else cout<<endl;
        cout<<mpi_rank<<":  unpack ranges:"<<endl;
        i_f = 0;
        for (range = prcomm->unpackRangeList.begin(); range!=prcomm->unpackRangeList.end(); range++)
        {
          cout<<mpi_rank<<":    -------------------------------- "<<endl;
          cout<<mpi_rank<<":    range->getBits()             : ";
          int ii;
          for (ii = 0; ii<32; ii++)
          {
            if (range->getBits()&(1<<ii))
            {
              cout<<"1";
            }
            else
            {
              cout<<"0";
            }
          }
          cout<<endl;
          cout<<mpi_rank<<":    range->getIndexFirst/Last    : "<<range->getIndexFirst()<<" "<<range->getIndexLast()
              <<endl;
          assert(range->getIndexFirst()==i_f);
          i_f = range->getIndexLast()+1;
        }
        assert(prcomm->unpackIndexVec.size()==i_f);
        cout<<mpi_rank<<":    -------------------------------- "<<endl;
        cout<<mpi_rank<<":  end of unpack ranges"<<endl;
      }
      cout<<mpi_rank<<":  =================================== "<<endl;
      cout<<mpi_rank<<":end of face prcomms"<<endl;

    }

    MPI_Pause("press return...");
  }

}

void Gp::dumpNodePrcommData()
{

  int rank;
  for (rank = 0; rank<mpi_size; rank++)
  {

    if (mpi_rank==rank)
    {

      cout<<mpi_rank<<":dumpNodePrcommData()..."<<endl;
      cout<<mpi_rank<<":nodePrcommList.size(): "<<nodePrcommList.size()<<endl;
      cout<<mpi_rank<<":node prcomms:"<<endl;
      list<Prcomm>::iterator prcomm;
      for (prcomm = nodePrcommList.begin(); prcomm!=nodePrcommList.end(); prcomm++)
      {
        cout<<mpi_rank<<":  =================================== "<<endl;
        cout<<mpi_rank<<":  prcomm->getNbrRank()           : "<<prcomm->getNbrRank()<<endl;
        cout<<mpi_rank<<":  prcomm->packIndexVec.size()    : "<<prcomm->packIndexVec.size()<<endl;
        cout<<mpi_rank<<":  prcomm->packIndexVec           : ";
        int i;
        for (i = 0; i<min(12, (int) prcomm->packIndexVec.size()); i++)
          cout<<prcomm->packIndexVec[i]<<" ";
        if (prcomm->packIndexVec.size()>12) cout<<mpi_rank<<":..."<<endl;
        else cout<<endl;
        cout<<mpi_rank<<":  prcomm->packRangeList.size()   : "<<prcomm->packRangeList.size()<<endl;
        cout<<mpi_rank<<":  pack ranges:"<<endl;
        int i_f = 0;
        list<Range>::iterator range;
        for (range = prcomm->packRangeList.begin(); range!=prcomm->packRangeList.end(); range++)
        {
          cout<<mpi_rank<<":    -------------------------------- "<<endl;
          cout<<mpi_rank<<":    range->getBits()             : "<<range->getBits()<<endl;
          cout<<mpi_rank<<":    range->getIndexFirst/Last    : "<<range->getIndexFirst()<<" "<<range->getIndexLast()
              <<endl;
          assert(range->getIndexFirst()==i_f);
          i_f = range->getIndexLast()+1;
        }
        assert(prcomm->packIndexVec.size()==i_f);
        cout<<mpi_rank<<":    -------------------------------- "<<endl;
        cout<<mpi_rank<<":  end of pack ranges"<<endl;
        cout<<mpi_rank<<":  prcomm->unpackRangeList.size() : "<<prcomm->unpackRangeList.size()<<endl;
        cout<<mpi_rank<<":  prcomm->unpackIndexVec.size()  : "<<prcomm->unpackIndexVec.size()<<endl;
        cout<<mpi_rank<<":  prcomm->unpackIndexVec         : ";
        for (i = 0; i<min(12, (int) prcomm->unpackIndexVec.size()); i++)
          cout<<prcomm->unpackIndexVec[i]<<" ";
        if (prcomm->unpackIndexVec.size()>12) cout<<mpi_rank<<":..."<<endl;
        else cout<<endl;
        cout<<mpi_rank<<":  unpack ranges:"<<endl;
        i_f = 0;
        for (range = prcomm->unpackRangeList.begin(); range!=prcomm->unpackRangeList.end(); range++)
        {
          cout<<mpi_rank<<":    -------------------------------- "<<endl;
          cout<<mpi_rank<<":    range->getBits()             : "<<range->getBits()<<endl;
          cout<<mpi_rank<<":    range->getIndexFirst/Last    : "<<range->getIndexFirst()<<" "<<range->getIndexLast()
              <<endl;
          assert(range->getIndexFirst()==i_f);
          i_f = range->getIndexLast()+1;
        }
        assert(prcomm->unpackIndexVec.size()==i_f);
        cout<<mpi_rank<<":    -------------------------------- "<<endl;
        cout<<mpi_rank<<":  end of unpack ranges"<<endl;
      }
      cout<<mpi_rank<<":  =================================== "<<endl;
      cout<<mpi_rank<<":end of node prcomms"<<endl;

    }

    MPI_Pause("press return...");
  }

}

void Gp::dump()
{

  cout<<"Gp::dump()"<<endl;

}

