#ifdef USENCCL
! Copyright (c) 2017, NVIDIA CORPORATION. All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions
! are met:
!  * Redistributions of source code must retain the above copyright
!    notice, this list of conditions and the following disclaimer.
!  * Redistributions in binary form must reproduce the above copyright
!    notice, this list of conditions and the following disclaimer in the
!    documentation and/or other materials provided with the distribution.
!  * Neither the name of NVIDIA CORPORATION nor the names of its
!    contributors may be used to endorse or promote products derived
!    from this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS ``AS IS'' AND ANY
! EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
! PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
! CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
! EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
! PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
! PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
! OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
! (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

!Start defines
#define NCCL_UNIQUE_ID_BYTES 128
!End defines

module nccl2For
use iso_c_binding
use cudaFor
implicit none
private
public :: ncclUniqueId
public :: ncclComm
public :: ncclResult,                 &
          ncclSuccess,                &
          ncclUnhandledCudaError,     &
          ncclSystemError,            &
          ncclInternalError,          &
          ncclInvalidDevicePointer,   &
          ncclInvalidRank,            &
          ncclUnsupportedDeviceCount, &
          ncclDeviceNotFound,         &
          ncclInvalidDeviceIndex,     &
          ncclLibWrapperNotSet,       &
          ncclCudaMallocFailed,       &
          ncclRankMismatch,           &
          ncclInvalidArgument,        &
          ncclInvalidType,            &
          ncclInvalidOperation,       &
          nccl_NUM_RESULTS
public :: ncclDataType, &
          ncclChar,     &
          ncclInt,      &
          ncclHalf,     &
          ncclFloat,    &
          ncclDouble,   &
          ncclInt64,    &
          ncclUInt64,   &
          nccl_NUM_TYPES
public :: ncclRedOp, &
          ncclSum,   &
          ncclProd,  &
          ncclMax,   &
          ncclMin,   &
          nccl_NUM_OPS
public :: ncclGetUniqueId
public :: ncclCommInitRank
public :: ncclGroupStart
public :: ncclGroupEnd

public :: ncclAllReduce
public :: ncclAllGather
public :: ncclReduce
public :: ncclBcast

#ifdef USENCCLP2P
public :: ncclSend
public :: ncclRecv
#endif

!Start types

!Start ncclUniqueId
type, bind(c) :: ncclUniqueId
character(c_char) :: internal(NCCL_UNIQUE_ID_BYTES)
end type ncclUniqueId
!End ncclUniqueId

!Start ncclComm
type, bind(c) :: ncclComm
type(c_ptr) :: member
end type ncclComm
!End ncclComm

!Start ncclResult
type, bind(c) :: ncclResult
integer(c_int) :: member
end type ncclResult

type(ncclResult), parameter :: ncclSuccess                = ncclResult( 0), &
                               ncclUnhandledCudaError     = ncclResult( 1), &
                               ncclSystemError            = ncclResult( 2), &
                               ncclInternalError          = ncclResult( 3), &
                               ncclInvalidDevicePointer   = ncclResult( 4), &
                               ncclInvalidRank            = ncclResult( 5), &
                               ncclUnsupportedDeviceCount = ncclResult( 6), &
                               ncclDeviceNotFound         = ncclResult( 7), &
                               ncclInvalidDeviceIndex     = ncclResult( 8), &
                               ncclLibWrapperNotSet       = ncclResult( 9), &
                               ncclCudaMallocFailed       = ncclResult(10), &
                               ncclRankMismatch           = ncclResult(11), &
                               ncclInvalidArgument        = ncclResult(12), &
                               ncclInvalidType            = ncclResult(13), &
                               ncclInvalidOperation       = ncclResult(14), &
                               nccl_NUM_RESULTS           = ncclResult(15)
!End ncclResult

!Start ncclDataType
type, bind(c) :: ncclDataType
integer(c_int) :: member
end type ncclDataType

type(ncclDataType), parameter :: ncclChar       = ncclDataType(0), &
                                 ncclInt        = ncclDataType(2), &
                                 ncclHalf       = ncclDataType(6), &
                                 ncclFloat      = ncclDataType(7), &
                                 ncclDouble     = ncclDataType(8), &
                                 ncclInt64      = ncclDataType(4), &
                                 ncclUInt64     = ncclDataType(5), &
                                 nccl_NUM_TYPES = ncclDataType(9)
!End ncclDataType

!Start ncclRedOp
type, bind(c) :: ncclRedOp
integer(c_int) :: member
end type ncclRedOp

type(ncclRedOp), parameter :: ncclSum      = ncclRedOp(0), &
                              ncclProd     = ncclRedOp(1), &
                              ncclMax      = ncclRedOp(2), &
                              ncclMin      = ncclRedOp(3), &
                              nccl_NUM_OPS = ncclRedOp(4)
!End ncclRedOp

!End types

!Start interfaces

!Start ncclGetUniqueId
interface ncclGetUniqueId
type(ncclResult) function ncclGetUniqueId(uniqueId) bind(c, name = 'ncclGetUniqueId')
import :: ncclResult, ncclUniqueId
implicit none
type(ncclUniqueId) :: uniqueId
end function ncclGetUniqueId
end interface ncclGetUniqueId
!End ncclGetUniqueId

!Start ncclCommInitRank
interface ncclCommInitRank
type(ncclResult) function ncclCommInitRank(comm, ndev, commId, rank) bind(c, name = 'ncclCommInitRank')
import :: c_int
import :: ncclResult, ncclUniqueId, ncclComm
implicit none
type(ncclComm) :: comm
integer(c_int), value :: ndev
type(ncclUniqueId), value :: commId
integer(c_int), value :: rank
end function ncclCommInitRank
end interface ncclCommInitRank
!End ncclCommInitRank

!Start ncclAllReduce
interface ncclAllReduce
type(ncclResult) function ncclAllReduce(sendbuff, recvbuff, count, &
                                        datatype, op, comm, stream) &
                          bind(c, name = 'ncclAllReduce')
import :: c_int
import :: c_devptr, cuda_stream_kind
import :: ncclResult, ncclComm, ncclDataType, ncclRedOp
implicit none
type(c_devptr), value :: sendbuff
type(c_devptr), value :: recvbuff
integer(c_int), value :: count
!integer(c_int), value :: datatype
!integer(c_int), value :: op
type(ncclDataType), value :: datatype
type(ncclRedOp), value :: op
type(ncclComm), value :: comm
integer(cuda_stream_kind), value :: stream
end function ncclAllReduce
end interface ncclAllReduce
!End ncclAllReduce

!Start ncclReduce
interface ncclReduce
type(ncclResult) function ncclReduce(sendbuff, recvbuff, count, &
                                        datatype, op, root, comm, stream) &
                          bind(c, name = 'ncclReduce')
import :: c_int
import :: c_devptr, cuda_stream_kind
import :: ncclResult, ncclComm, ncclDataType, ncclRedOp
implicit none
type(c_devptr), value :: sendbuff
type(c_devptr), value :: recvbuff
integer(c_int), value :: count
!integer(c_int), value :: datatype
!integer(c_int), value :: op
type(ncclDataType), value :: datatype
type(ncclRedOp), value :: op
integer(c_int), value :: root
type(ncclComm), value :: comm
integer(cuda_stream_kind), value :: stream
end function ncclReduce
end interface ncclReduce
!End ncclReduce

!Start ncclBcast
interface ncclBcast
type(ncclResult) function ncclBcast(buff, count, &
                                        datatype, root, comm, stream) &
                          bind(c, name = 'ncclBcast')
import :: c_int
import :: c_devptr, cuda_stream_kind
import :: ncclResult, ncclComm, ncclDataType, ncclRedOp
implicit none
type(c_devptr), value :: buff
integer(c_int), value :: count
!integer(c_int), value :: datatype
!integer(c_int), value :: op
type(ncclDataType), value :: datatype
integer(c_int), value :: root
type(ncclComm), value :: comm
integer(cuda_stream_kind), value :: stream
end function ncclBcast
end interface ncclBcast
!End ncclBcast

!Start ncclAllGather
interface ncclAllGather
type(ncclResult) function ncclAllGather(sendbuff, recvbuff, count, &
                                        datatype, comm, stream) &
                          bind(c, name = 'ncclAllGather')
import :: c_int
import :: c_devptr, cuda_stream_kind
import :: ncclResult, ncclComm, ncclDataType
implicit none
type(c_devptr), value :: sendbuff
type(c_devptr), value :: recvbuff
integer(c_int), value :: count
type(ncclDataType), value :: datatype
type(ncclComm), value :: comm
integer(cuda_stream_kind), value :: stream
end function ncclAllGather
end interface ncclAllGather
!End ncclReduce

#ifdef USENCCLP2P
!Start ncclSend
interface ncclSend
type(ncclResult) function ncclSend(sendbuff, count, datatype, peer, comm, stream) &
                          bind(c, name = 'ncclSend')
import :: c_int
import :: c_devptr, cuda_stream_kind
import :: ncclResult, ncclComm, ncclDataType
implicit none
type(c_devptr), value :: sendbuff
integer(c_int), value :: count
integer(c_int), value :: peer
type(ncclDataType), value :: datatype
type(ncclComm), value :: comm
integer(cuda_stream_kind), value :: stream
end function ncclSend
end interface ncclSend
!End ncclSend

!Start ncclRecv
interface ncclRecv
type(ncclResult) function ncclRecv(recvbuff, count, datatype, peer, comm, stream) &
                          bind(c, name = 'ncclRecv')
import :: c_int
import :: c_devptr, cuda_stream_kind
import :: ncclResult, ncclComm, ncclDataType
implicit none
type(c_devptr), value :: recvbuff
integer(c_int), value :: count
integer(c_int), value :: peer
type(ncclDataType), value :: datatype
type(ncclComm), value :: comm
integer(cuda_stream_kind), value :: stream
end function ncclRecv
end interface ncclRecv
!End ncclRecv
#endif


!Start ncclGroupStart
interface ncclGroupStart
type(ncclResult) function ncclGroupStart() bind(c, name = 'ncclGroupStart')
end function ncclGroupStart
end interface ncclGroupStart
!End ncclGroupStart

!End ncclGroupEnd
interface ncclGroupEnd
type(ncclResult) function ncclGroupEnd() bind(c, name = 'ncclGroupEnd')
end function ncclGroupEnd
end interface ncclGroupEnd
!End ncclGroupEnd

end module nccl2For
!End nccl module
#else
module nccl2For
contains
  subroutine nccl2For_dummy
    write(*,*)'Im a DEC compiler so I need this line'
  end subroutine nccl2For_dummy
end module nccl2For
#endif
