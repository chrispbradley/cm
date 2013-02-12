!> \file
!> \author Chris Bradley
!> \brief Implements lists of base types.
!> \todo Fix up and have this module use the sorting module for sorts
!>
!> \section LICENSE
!>
!> Version: MPL 1.1/GPL 2.0/LGPL 2.1
!>
!> The contents of this file are subject to the Mozilla Public License
!> Version 1.1 (the "License"); you may not use this file except in
!> compliance with the License. You may obtain a copy of the License at
!> http://www.mozilla.org/MPL/
!>
!> Software distributed under the License is distributed on an "AS IS"
!> basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See the
!> License for the specific language governing rights and limitations
!> under the License.
!>
!> The Original Code is OpenCMISS
!>
!> The Initial Developer of the Original Code is University of Auckland,
!> Auckland, New Zealand, the University of Oxford, Oxford, United
!> Kingdom and King's College, London, United Kingdom. Portions created
!> by the University of Auckland, the University of Oxford and King's
!> College, London are Copyright (C) 2007-2010 by the University of
!> Auckland, the University of Oxford and King's College, London.
!> All Rights Reserved.
!>
!> Contributor(s):
!>
!> Alternatively, the contents of this file may be used under the terms of
!> either the GNU General Public License Version 2 or later (the "GPL"), or
!> the GNU Lesser General Public License Version 2.1 or later (the "LGPL"),
!> in which case the provisions of the GPL or the LGPL are applicable instead
!> of those above. If you wish to allow use of your version of this file only
!> under the terms of either the GPL or the LGPL, and not to allow others to
!> use your version of this file under the terms of the MPL, indicate your
!> decision by deleting the provisions above and replace them with the notice
!> and other provisions required by the GPL or the LGPL. If you do not delete
!> the provisions above, a recipient may use your version of this file under
!> the terms of any one of the MPL, the GPL or the LGPL.
!>

!> Implements lists of base types.
MODULE Lists

  USE BASE_ROUTINES
  USE CONSTANTS
  USE ISO_VARYING_STRING
  USE KINDS
  USE STRINGS
  USE TYPES
  
  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !> \addtogroup LISTS_DataType LISTS::DataType
  !> \brief Data type parameters for a list.
  !> \see LISTS
  !>@{
  INTEGER(INTG), PARAMETER :: LIST_INTG_TYPE=INTEGER_TYPE !<Integer data type for a list \see LISTS_DataType,LISTS
  INTEGER(INTG), PARAMETER :: LIST_SP_TYPE=SINGLE_REAL_TYPE !<Single precision real data type for a list \see LISTS_DataType,LISTS
  INTEGER(INTG), PARAMETER :: LIST_DP_TYPE=DOUBLE_REAL_TYPE !<Double precision real data type for a list \see LISTS_DataType,LISTS
  !>@}
  
  !> \addtogroup LISTS_SortingOrder LISTS::SortingOrder
  !> \brief Sorting order parameters for a list.
  !> \see LISTS
  !>@{
  INTEGER(INTG), PARAMETER :: LIST_UNSORTED_TYPE=1 !<Unsorted list type \see LISTS_SortingOrder,LISTS
  INTEGER(INTG), PARAMETER :: LIST_SORT_ASCENDING_TYPE=2 !<Ascending order for sort \see LISTS_SortingOrder,LISTS
  INTEGER(INTG), PARAMETER :: LIST_SORT_DESCENDING_TYPE=3 !<Descending order for sort \see LISTS_SortingOrder,LISTS
  !>@}

  !> \addtogroup LISTSSortingMethod LISTS::SortingMethod
  !> \brief Sorting method parameters for a list.
  !> \see LISTS
  !>@{
  INTEGER(INTG), PARAMETER :: LIST_BUBBLE_SORT_METHOD=1 !<Bubble sort method \see LISTS_SortingMethod,LISTS
  INTEGER(INTG), PARAMETER :: LIST_SHELL_SORT_METHOD=2 !<Shell sort method \see LISTS_SortingMethod,LISTS
  INTEGER(INTG), PARAMETER :: LIST_HEAP_SORT_METHOD=3 !<Heap sort method \see LISTS_SortingMethod,LISTS
  !>@}

  !Module types

 !Module variables
  
  !Interfaces

  !>Adds an item to the end of a list \see LISTS.
  INTERFACE List_ItemAdd
    MODULE PROCEDURE LIST_ITEM_ADD_INTG1
    MODULE PROCEDURE LIST_ITEM_ADD_INTG2
    MODULE PROCEDURE LIST_ITEM_ADD_SP1
    MODULE PROCEDURE LIST_ITEM_ADD_SP2
    MODULE PROCEDURE LIST_ITEM_ADD_DP1
    MODULE PROCEDURE LIST_ITEM_ADD_DP2
  END INTERFACE List_ItemAdd
  
  !>Adds an item to the end of a list \see LISTS.
  INTERFACE LIST_ITEM_ADD
    MODULE PROCEDURE LIST_ITEM_ADD_INTG1
    MODULE PROCEDURE LIST_ITEM_ADD_INTG2
    MODULE PROCEDURE LIST_ITEM_ADD_SP1
    MODULE PROCEDURE LIST_ITEM_ADD_SP2
    MODULE PROCEDURE LIST_ITEM_ADD_DP1
    MODULE PROCEDURE LIST_ITEM_ADD_DP2
  END INTERFACE !LIST_ITEM_ADD
  
  !>Sets an item in the list \see LISTS.
  INTERFACE List_ItemSet
    MODULE PROCEDURE LIST_ITEM_SET_INTG1
    MODULE PROCEDURE LIST_ITEM_SET_INTG2
    MODULE PROCEDURE LIST_ITEM_SET_SP1
    MODULE PROCEDURE LIST_ITEM_SET_SP2
    MODULE PROCEDURE LIST_ITEM_SET_DP1
    MODULE PROCEDURE LIST_ITEM_SET_DP2
  END INTERFACE List_ItemSet
  
  !>Sets an item in the list \see LISTS.
  INTERFACE LIST_ITEM_SET
    MODULE PROCEDURE LIST_ITEM_SET_INTG1
    MODULE PROCEDURE LIST_ITEM_SET_INTG2
    MODULE PROCEDURE LIST_ITEM_SET_SP1
    MODULE PROCEDURE LIST_ITEM_SET_SP2
    MODULE PROCEDURE LIST_ITEM_SET_DP1
    MODULE PROCEDURE LIST_ITEM_SET_DP2
  END INTERFACE !LIST_ITEM_SET
  
  !>Returns an item in a list at a specififed position. \see LISTS.
  INTERFACE List_ItemGet
    MODULE PROCEDURE LIST_ITEM_GET_INTG1
    MODULE PROCEDURE LIST_ITEM_GET_INTG2
    MODULE PROCEDURE LIST_ITEM_GET_SP1
    MODULE PROCEDURE LIST_ITEM_GET_SP2
    MODULE PROCEDURE LIST_ITEM_GET_DP1
    MODULE PROCEDURE LIST_ITEM_GET_DP2
  END INTERFACE List_ItemGet

  !>Returns an item in a list at a specififed position. \see LISTS.
  INTERFACE LIST_ITEM_GET
    MODULE PROCEDURE LIST_ITEM_GET_INTG1
    MODULE PROCEDURE LIST_ITEM_GET_INTG2
    MODULE PROCEDURE LIST_ITEM_GET_SP1
    MODULE PROCEDURE LIST_ITEM_GET_SP2
    MODULE PROCEDURE LIST_ITEM_GET_DP1
    MODULE PROCEDURE LIST_ITEM_GET_DP2
  END INTERFACE !LIST_ITEM_GET

  !>Determines if an item is in a list and returns the position of the item \see LISTS.
  INTERFACE List_ItemInList
    MODULE PROCEDURE LIST_ITEM_IN_LIST_INTG1
    MODULE PROCEDURE LIST_ITEM_IN_LIST_SP1
    MODULE PROCEDURE LIST_ITEM_IN_LIST_DP1
  END INTERFACE List_ItemInList

  !>Determines if an item is in a list and returns the position of the item \see LISTS.
  INTERFACE LIST_ITEM_IN_LIST
    MODULE PROCEDURE LIST_ITEM_IN_LIST_INTG1
    MODULE PROCEDURE LIST_ITEM_IN_LIST_SP1
    MODULE PROCEDURE LIST_ITEM_IN_LIST_DP1
  END INTERFACE !LIST_ITEM_IN_LIST

  !>Detaches the list values from a list and returns them as a pointer to a array of base type before destroying the list \see LISTS.
  INTERFACE List_DetachAndDestroy
    MODULE PROCEDURE LIST_DETACH_AND_DESTROY_INTG1
    MODULE PROCEDURE LIST_DETACH_AND_DESTROY_INTG2
    MODULE PROCEDURE LIST_DETACH_AND_DESTROY_SP1
    MODULE PROCEDURE LIST_DETACH_AND_DESTROY_SP2
    MODULE PROCEDURE LIST_DETACH_AND_DESTROY_DP1
    MODULE PROCEDURE LIST_DETACH_AND_DESTROY_DP2
  END INTERFACE List_DetachAndDestroy

  !>Detaches the list values from a list and returns them as a pointer to a array of base type before destroying the list \see LISTS.
  INTERFACE LIST_DETACH_AND_DESTROY
    MODULE PROCEDURE LIST_DETACH_AND_DESTROY_INTG1
    MODULE PROCEDURE LIST_DETACH_AND_DESTROY_INTG2
    MODULE PROCEDURE LIST_DETACH_AND_DESTROY_SP1
    MODULE PROCEDURE LIST_DETACH_AND_DESTROY_SP2
    MODULE PROCEDURE LIST_DETACH_AND_DESTROY_DP1
    MODULE PROCEDURE LIST_DETACH_AND_DESTROY_DP2
  END INTERFACE !LIST_DETACH_AND_DESTROY

  !>Searches a list for a given value and returns the position in the list if the value exists \see LISTS.
  INTERFACE LIST_SEARCH
    MODULE PROCEDURE LIST_SEARCH_INTG_ARRAY
    MODULE PROCEDURE LIST_SEARCH_SP_ARRAY
    MODULE PROCEDURE LIST_SEARCH_DP_ARRAY
  END INTERFACE !LIST_SEARCH

  !>Searches a list using the linear search method.
  INTERFACE List_SearchLinear
    MODULE PROCEDURE LIST_SEARCH_LINEAR_INTG_ARRAY
    MODULE PROCEDURE LIST_SEARCH_LINEAR_SP_ARRAY
    MODULE PROCEDURE LIST_SEARCH_LINEAR_DP_ARRAY
  END INTERFACE List_SearchLinear

  !>Searches a list using the linear search method.
  INTERFACE LIST_SEARCH_LINEAR
    MODULE PROCEDURE LIST_SEARCH_LINEAR_INTG_ARRAY
    MODULE PROCEDURE LIST_SEARCH_LINEAR_SP_ARRAY
    MODULE PROCEDURE LIST_SEARCH_LINEAR_DP_ARRAY
  END INTERFACE !LIST_SEARCH_LINEAR

  !>Sorts a list into ascending order.
  INTERFACE LIST_SORT
    MODULE PROCEDURE LIST_SORT_LIST
    MODULE PROCEDURE LIST_SORT_INTG1_ARRAY
    MODULE PROCEDURE LIST_SORT_INTG2_ARRAY
    MODULE PROCEDURE LIST_SORT_SP1_ARRAY
    MODULE PROCEDURE LIST_SORT_SP2_ARRAY
    MODULE PROCEDURE LIST_SORT_DP1_ARRAY
    MODULE PROCEDURE LIST_SORT_DP2_ARRAY
  END INTERFACE !LIST_SORT

  !>Sorts a list into assending order using the bubble sort method.
  INTERFACE List_SortBubble
    MODULE PROCEDURE LIST_SORT_BUBBLE_INTG1_ARRAY
    MODULE PROCEDURE LIST_SORT_BUBBLE_INTG2_ARRAY
    MODULE PROCEDURE LIST_SORT_BUBBLE_SP1_ARRAY
    MODULE PROCEDURE LIST_SORT_BUBBLE_SP2_ARRAY
    MODULE PROCEDURE LIST_SORT_BUBBLE_DP1_ARRAY
    MODULE PROCEDURE LIST_SORT_BUBBLE_DP2_ARRAY
  END INTERFACE List_SortBubble

  !>Sorts a list into assending order using the bubble sort method.
  INTERFACE LIST_SORT_BUBBLE
    MODULE PROCEDURE LIST_SORT_BUBBLE_INTG1_ARRAY
    MODULE PROCEDURE LIST_SORT_BUBBLE_INTG2_ARRAY
    MODULE PROCEDURE LIST_SORT_BUBBLE_SP1_ARRAY
    MODULE PROCEDURE LIST_SORT_BUBBLE_SP2_ARRAY
    MODULE PROCEDURE LIST_SORT_BUBBLE_DP1_ARRAY
    MODULE PROCEDURE LIST_SORT_BUBBLE_DP2_ARRAY
  END INTERFACE !LIST_SORT_BUBBLE

  !>Sorts a list into assending order using the heap sort method.
  INTERFACE List_SortHeap
    MODULE PROCEDURE LIST_SORT_HEAP_INTG1_ARRAY
    MODULE PROCEDURE LIST_SORT_HEAP_INTG2_ARRAY
    MODULE PROCEDURE LIST_SORT_HEAP_SP1_ARRAY
    MODULE PROCEDURE LIST_SORT_HEAP_SP2_ARRAY
    MODULE PROCEDURE LIST_SORT_HEAP_DP1_ARRAY
    MODULE PROCEDURE LIST_SORT_HEAP_DP2_ARRAY
  END INTERFACE List_SortHeap

  !>Sorts a list into assending order using the heap sort method.
  INTERFACE LIST_SORT_HEAP
    MODULE PROCEDURE LIST_SORT_HEAP_INTG1_ARRAY
    MODULE PROCEDURE LIST_SORT_HEAP_INTG2_ARRAY
    MODULE PROCEDURE LIST_SORT_HEAP_SP1_ARRAY
    MODULE PROCEDURE LIST_SORT_HEAP_SP2_ARRAY
    MODULE PROCEDURE LIST_SORT_HEAP_DP1_ARRAY
    MODULE PROCEDURE LIST_SORT_HEAP_DP2_ARRAY
  END INTERFACE !LIST_SORT_HEAP

  !>Sorts a list into either assending or descending order using the shell sort method.
  INTERFACE List_SortShell
    MODULE PROCEDURE LIST_SORT_SHELL_INTG1_ARRAY
    MODULE PROCEDURE LIST_SORT_SHELL_INTG2_ARRAY
    MODULE PROCEDURE LIST_SORT_SHELL_SP1_ARRAY
    MODULE PROCEDURE LIST_SORT_SHELL_SP2_ARRAY
    MODULE PROCEDURE LIST_SORT_SHELL_DP1_ARRAY
    MODULE PROCEDURE LIST_SORT_SHELL_DP2_ARRAY
  END INTERFACE List_SortShell

  !>Sorts a list into either assending or descending order using the shell sort method.
  INTERFACE LIST_SORT_SHELL
    MODULE PROCEDURE LIST_SORT_SHELL_INTG1_ARRAY
    MODULE PROCEDURE LIST_SORT_SHELL_INTG2_ARRAY
    MODULE PROCEDURE LIST_SORT_SHELL_SP1_ARRAY
    MODULE PROCEDURE LIST_SORT_SHELL_SP2_ARRAY
    MODULE PROCEDURE LIST_SORT_SHELL_DP1_ARRAY
    MODULE PROCEDURE LIST_SORT_SHELL_DP2_ARRAY
  END INTERFACE !LIST_SORT_SHELL

  !>Calculates the intersection of two arrays
  INTERFACE LIST_INTERSECTION
    MODULE PROCEDURE LIST_INTERSECTION_INTG_ARRAY
  END INTERFACE LIST_INTERSECTION

  !>Checks whether an array is a subset of another array
  INTERFACE List_SubsetOf
    MODULE PROCEDURE LISTS_SUBSET_OF_INTG_ARRAY
  END INTERFACE List_SubsetOf

 !>Checks whether an array is a subset of another array
  INTERFACE LIST_SUBSET_OF
    MODULE PROCEDURE LISTS_SUBSET_OF_INTG_ARRAY
  END INTERFACE LIST_SUBSET_OF

  PUBLIC LIST_INTG_TYPE,LIST_SP_TYPE,LIST_DP_TYPE

  PUBLIC List_AppendList

  PUBLIC List_ClearItems

  PUBLIC List_CreateFinish,List_CreateStart

  PUBLIC LIST_CREATE_FINISH,LIST_CREATE_START

  PUBLIC List_MutableSet
  
  PUBLIC LIST_MUTABLE_SET

  PUBLIC List_DataDimensionSet

  PUBLIC LIST_DATA_DIMENSION_SET

  PUBLIC List_DataTypeSet
  
  PUBLIC LIST_DATA_TYPE_SET

  PUBLIC List_Destroy

  PUBLIC List_DetachAndDestroy

  PUBLIC LIST_DETACH_AND_DESTROY

  PUBLIC List_InitialSizeSet

  PUBLIC LIST_INITIAL_SIZE_SET

  PUBLIC List_ItemAdd
  
  PUBLIC LIST_ITEM_ADD

  PUBLIC List_ItemSet
  
  PUBLIC LIST_ITEM_SET

  PUBLIC List_ItemDelete
  
  PUBLIC LIST_ITEM_DELETE

  PUBLIC List_ItemGet

  PUBLIC LIST_ITEM_GET

  PUBLIC List_KeyDimensionSet

  PUBLIC List_NumberOfItemsGet

  PUBLIC LIST_NUMBER_OF_ITEMS_GET

  PUBLIC List_RemoveDuplicates

  PUBLIC LIST_REMOVE_DUPLICATES

  PUBLIC List_Search

  PUBLIC List_SearchLinear

  PUBLIC LIST_SEARCH_LINEAR

  PUBLIC List_Sort

  PUBLIC List_SortBubble,List_SortHeap,List_SortShell
  
  PUBLIC LIST_SORT_BUBBLE,LIST_SORT_HEAP,LIST_SORT_SHELL

  PUBLIC List_Intersection

  PUBLIC List_SubsetOf

  PUBLIC LIST_SUBSET_OF

  PUBLIC List_ItemInList
  
  PUBLIC LIST_ITEM_IN_LIST

CONTAINS

  !
  !================================================================================================================================
  !

  !>Finishes the creation of a list created with List_CreateStart \see Lists::List_CreateStart
  SUBROUTINE List_CreateFinish(list,err,error,*)

    !Argument Variables
    TYPE(LIST_TYPE), POINTER, INTENT(INOUT) :: list !<A pointer to the list to finish
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError,localError
    
    CALL Enters("List_CreateFinish",err,error,*998)

    IF(ASSOCIATED(list)) THEN
      IF(list%LIST_FINISHED) THEN
        CALL FlagError("List is already finished.",err,error,*998)
      ELSE
        !Allocate the list
        IF(list%DATA_DIMENSION==1) THEN
          SELECT CASE(list%DATA_TYPE)
          CASE(LIST_INTG_TYPE)
            ALLOCATE(list%LIST_INTG(list%INITIAL_SIZE),STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate list integer data.",err,error,*999)
          CASE(LIST_SP_TYPE)
            ALLOCATE(list%LIST_SP(list%INITIAL_SIZE),STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate list single precision data.",err,error,*999)
          CASE(LIST_DP_TYPE)
            ALLOCATE(list%LIST_DP(list%INITIAL_SIZE),STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate list double precision data.",err,error,*999)
          CASE DEFAULT
            localError="The list data type of "//TRIM(NumberToVString(list%DATA_TYPE,"*",err,error))//" is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        ELSE
          SELECT CASE(list%DATA_TYPE)
          CASE(LIST_INTG_TYPE)
            ALLOCATE(list%LIST_INTG2(list%DATA_DIMENSION,list%INITIAL_SIZE),STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate list integer data.",err,error,*999)
          CASE(LIST_SP_TYPE)
            ALLOCATE(list%LIST_SP2(list%DATA_DIMENSION,list%INITIAL_SIZE),STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate list single precision data.",err,error,*999)
          CASE(LIST_DP_TYPE)
            ALLOCATE(list%LIST_DP2(list%DATA_DIMENSION,list%INITIAL_SIZE),STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate list double precision data.",err,error,*999)
          CASE DEFAULT
           localError="The list data type of "//TRIM(NumberToVString(list%DATA_TYPE,"*",err,error))//" is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        ENDIF
        list%size=list%INITIAL_SIZE
        list%LIST_FINISHED=.TRUE.
      ENDIF
    ELSE
      CALL FlagError("List is not associated.",err,error,*998)
    ENDIF

    Call Exits("List_CreateFinish")
    RETURN
999 CALL List_Finalise(list,dummyErr,dummyError,*998)
998 CALL Errors("List_CreateFinish",err,error)
    Call Exits("List_CreateFinish")
    RETURN 1
  END SUBROUTINE List_CreateFinish

  !
  !================================================================================================================================
  !

  !>Finishes the creation of a list created with LIST_CREATE_START \see LISTS::LIST_CREATE_START.
  SUBROUTINE LIST_CREATE_FINISH(LIST,err,error,*)

    !Argument Variables
    TYPE(LIST_TYPE), POINTER, INTENT(INOUT) :: LIST !<A pointer to the list to finish
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    !Local Variables

    CALL List_CreateFinish(LIST,err,error,*999)

    RETURN
999 RETURN 1    
  END SUBROUTINE LIST_CREATE_FINISH

  !
  !================================================================================================================================
  !

  !>Starts the creation of a list and returns a pointer to the created list \see Lists::List_CreateFinish.
  SUBROUTINE List_CreateStart(list,err,error,*)

    !Argument Variables
    TYPE(LIST_TYPE), POINTER, INTENT(INOUT) :: list !<On exit, pointer to the list to create. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code..
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    !Local Variables

    CALL Enters("List_CreateStart",err,error,*999)

    CALL List_Initialise(list,err,error,*999)
    
    Call Exits("List_CreateStart")
    RETURN
999 CALL Errors("List_CreateStart",err,error)
    CALL Exits("List_CreateStart")
    RETURN 1
  END SUBROUTINE List_CreateStart

  !
  !================================================================================================================================
  !

  !>Starts the creation of a list and returns a pointer to the created list \see LISTS::LIST_CREATE_FINISH.
  SUBROUTINE LIST_CREATE_START(LIST,err,error,*)

    !Argument Variables
    TYPE(LIST_TYPE), POINTER, INTENT(INOUT) :: LIST !<On exit, pointer to the list to create. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code..
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    !Local Variables

    CALL List_CreateStart(LIST,err,error,*999)
    
    RETURN
999 RETURN 1
  END SUBROUTINE LIST_CREATE_START

  !
  !================================================================================================================================
  !

  !>Sets/changes the data dimension for a list.
  SUBROUTINE List_DataDimensionSet(list,dataDimension,err,error,*)

    !Argument Variables
    TYPE(LIST_TYPE), POINTER, INTENT(IN) :: list !<A pointer to the list 
    INTEGER(INTG), INTENT(IN) :: dataDimension !<The data dimension of the list to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    CALL Enters("List_DataDimensionSet",err,error,*999)

    IF(ASSOCIATED(list)) THEN
      IF(list%LIST_FINISHED) THEN
        CALL FlagError("List has been finished.",err,error,*999)
      ELSE
        IF(dataDimension>0) THEN
          list%DATA_DIMENSION=dataDimension
        ELSE
          localError="The specified data dimension of "//TRIM(NumberToVString(dataDimension,"*",err,error))// &
            & " is invalid. The dimension must be > 0."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("List is not associated.",err,error,*999)
    ENDIF

    Call Exits("List_DataDimensionSet")
    RETURN
999 CALL Errors("List_DataDimensionSet",err,error)
    Call Exits("List_DataDimensionSet")
    RETURN 1
  END SUBROUTINE List_DataDimensionSet

  !
  !================================================================================================================================
  !

  !>Sets/changes the data dimension for a list.
  SUBROUTINE LIST_DATA_DIMENSION_SET(LIST,DATA_DIMENSION,err,error,*)

    !Argument Variables
    TYPE(LIST_TYPE), POINTER, INTENT(IN) :: LIST !<A pointer to the list 
    INTEGER(INTG), INTENT(IN) :: DATA_DIMENSION !<The data dimension of the list to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL List_DataDimensionSet(LIST,DATA_DIMENSION,err,error,*999)

    RETURN
999 RETURN 1
  END SUBROUTINE LIST_DATA_DIMENSION_SET

  !
  !================================================================================================================================
  !

  !>Sets/changes the data dimension for a list.
  SUBROUTINE List_MutableSet(list,mutable,err,error,*)

    !Argument Variables
    TYPE(LIST_TYPE), POINTER, INTENT(IN) :: list !<A pointer to the list 
    LOGICAL, INTENT(IN) :: mutable !<The mutability of the list to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    CALL Enters("List_MutableSet",err,error,*999)

    IF(ASSOCIATED(list)) THEN
      IF(list%LIST_FINISHED) THEN
        CALL FlagError("List has been finished.",err,error,*999)
      ELSE
        list%MUTABLE = mutable
      ENDIF
    ELSE
      CALL FlagError("List is not associated.",err,error,*999)
    ENDIF

    Call Exits("List_MutableSet")
    RETURN
999 CALL Errors("List_MutableSet",err,error)
    Call Exits("List_MutableSet")
    RETURN 1
  END SUBROUTINE List_MutableSet

  !
  !================================================================================================================================
  !

  !>Sets/changes the data dimension for a list.
  SUBROUTINE LIST_MUTABLE_SET(LIST,MUTABLE,err,error,*)

    !Argument Variables
    TYPE(LIST_TYPE), POINTER, INTENT(IN) :: LIST !<A pointer to the list 
    LOGICAL, INTENT(IN) :: MUTABLE !<The mutability of the list to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    CALL List_MutableSet(LIST,MUTABLE,err,error,*999)
    
    RETURN
999 RETURN 1
  END SUBROUTINE LIST_MUTABLE_SET

  !
  !================================================================================================================================
  !

  !>Sets/changes the data type for a list.
  SUBROUTINE List_DataTypeSet(list,dataType,err,error,*)

    !Argument Variables
    TYPE(LIST_TYPE), POINTER, INTENT(IN) :: list !<A pointer to the list 
    INTEGER(INTG), INTENT(IN) :: dataType !<The data type of the list to set \see Lists_DataType,Lists
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    CALL Enters("List_DataTypeSet",err,error,*999)

    IF(ASSOCIATED(list)) THEN
      IF(list%LIST_FINISHED) THEN
        CALL FlagError("List has been finished.",err,error,*999)
      ELSE
        SELECT CASE(dataType)
        CASE(LIST_INTG_TYPE)
          list%DATA_TYPE=LIST_INTG_TYPE
        CASE(LIST_SP_TYPE)
          list%DATA_TYPE=LIST_SP_TYPE
        CASE(LIST_DP_TYPE)
          list%DATA_TYPE=LIST_DP_TYPE
         CASE DEFAULT
          localError="The data type of "//TRIM(NumberToVString(dataType,"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ENDIF
    ELSE
      CALL FlagError("List is not associated.",err,error,*999)
    ENDIF

    Call Exits("List_DataTypeSet")
    RETURN
999 CALL Errors("List_DataTypeSet",err,error)
    Call Exits("List_DataTypeSet")
    RETURN 1
  END SUBROUTINE List_DataTypeSet

  !
  !================================================================================================================================
  !

  !>Sets/changes the data type for a list.
  SUBROUTINE LIST_DATA_TYPE_SET(LIST,DATA_TYPE,err,error,*)

    !Argument Variables
    TYPE(LIST_TYPE), POINTER, INTENT(IN) :: LIST !<A pointer to the list 
    INTEGER(INTG), INTENT(IN) :: DATA_TYPE !<The data type of the list to set \see LISTS_DataType,LISTS
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL List_DataTypeSet(LIST,DATA_TYPE,err,error,*999)
    
    RETURN
999 RETURN 1
  END SUBROUTINE LIST_DATA_TYPE_SET

  !
  !================================================================================================================================
  !

  !>Destroys a list.
  SUBROUTINE List_Destroy(list,err,error,*)

    !Argument Variables
    TYPE(LIST_TYPE), POINTER, INTENT(INOUT) :: list !<A pointer to the list to destroy
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    !Local Variables

    CALL Enters("List_Destroy",err,error,*999)

    IF(ASSOCIATED(list)) THEN
      CALL List_Finalise(list,err,error,*999)
    ELSE
      CALL FlagError("List is not associated.",err,error,*999)
    ENDIF

    Call Exits("List_Destroy")
    RETURN
999 CALL Errors("List_Destroy",err,error)
    Call Exits("List_Destroy")
    RETURN 1
  END SUBROUTINE List_Destroy

  !
  !================================================================================================================================
  !

  !>Finalises a list and deallocates all memory.
  SUBROUTINE List_Finalise(list,err,error,*)    

    !Argument Variables
    TYPE(LIST_TYPE), POINTER, INTENT(INOUT) :: list !<A pointer to the list to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    !Local Variables

    CALL Enters("List_Finalise",err,error,*999)

    IF(ASSOCIATED(list)) THEN
      IF(ALLOCATED(list%LIST_INTG)) DEALLOCATE(list%LIST_INTG)
      IF(ALLOCATED(list%LIST_INTG2)) DEALLOCATE(list%LIST_INTG2)
      IF(ALLOCATED(list%LIST_SP)) DEALLOCATE(list%LIST_SP)
      IF(ALLOCATED(list%LIST_SP2)) DEALLOCATE(list%LIST_SP2)
      IF(ALLOCATED(list%LIST_DP)) DEALLOCATE(list%LIST_DP)
      IF(ALLOCATED(list%LIST_DP2)) DEALLOCATE(list%LIST_DP2)
      DEALLOCATE(list)
    ENDIF

    Call Exits("List_Finalise")
    RETURN
999 CALL Errors("List_Finalise",err,error)
    Call Exits("List_Finalise")
    RETURN 1
  END SUBROUTINE List_Finalise

  !
  !================================================================================================================================
  !

  !>Appends a list to the end of this list
  SUBROUTINE List_AppendList(list,appendedList,err,error,*)
   !Argument Variables
    TYPE(LIST_TYPE), POINTER, INTENT(INOUT) :: list !<A pointer to the list
    TYPE(LIST_TYPE), POINTER, INTENT(IN) :: appendedList !<The list to append
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: newSize
    INTEGER(INTG), ALLOCATABLE :: newListIntg(:)
    REAL(SP), ALLOCATABLE :: newListSP(:)
    REAL(DP), ALLOCATABLE :: newListDP(:)
    INTEGER(C_INT), ALLOCATABLE :: newListCInt(:)
    TYPE(VARYING_STRING) :: localError

    CALL ENTERS("List_AppendList",err,error,*999)

    IF(ASSOCIATED(list)) THEN
      IF(list%LIST_FINISHED) THEN
        IF(ASSOCIATED(appendedList)) THEN
          IF(appendedList%LIST_FINISHED) THEN
            IF(list%DATA_TYPE==appendedList%DATA_TYPE) THEN
              IF(list%DATA_DIMENSION==appendedList%DATA_DIMENSION) THEN
                SELECT CASE(list%DATA_DIMENSION)
                CASE(1)
                  SELECT CASE(list%DATA_TYPE)
                  CASE(LIST_INTG_TYPE)
                    IF(list%SIZE<list%NUMBER_IN_LIST+appendedList%NUMBER_IN_LIST) THEN
                      !Reallocate
                      newSize=MAX(2*list%NUMBER_IN_LIST,list%NUMBER_IN_LIST+appendedList%NUMBER_IN_LIST*2)
                      ALLOCATE(newListIntg(newSize),stat=err)
                      IF(err/=0) CALL FLAG_ERROR("Could not allocate new list.",err,ERROR,*999)
                      newListIntg(1:list%NUMBER_IN_LIST)=list%LIST_INTG(1:list%NUMBER_IN_LIST)
                      CALL MOVE_ALLOC(newListIntg,list%LIST_INTG)
                      list%SIZE=newSize
                    END IF
                    list%LIST_INTG(list%NUMBER_IN_LIST+1:list%NUMBER_IN_LIST+appendedList%NUMBER_IN_LIST)= &
                      & appendedList%LIST_INTG(1:appendedList%NUMBER_IN_LIST)
                    list%NUMBER_IN_LIST=list%NUMBER_IN_LIST+appendedList%NUMBER_IN_LIST
                  CASE(LIST_SP_TYPE)
                    IF(list%SIZE<list%NUMBER_IN_LIST+appendedList%NUMBER_IN_LIST) THEN
                      !Reallocate
                      newSize=MAX(2*list%NUMBER_IN_LIST,list%NUMBER_IN_LIST+appendedList%NUMBER_IN_LIST*2)
                      ALLOCATE(newListSP(newSize),stat=err)
                      IF(err/=0) CALL FLAG_ERROR("Could not allocate new list.",err,ERROR,*999)
                      newListSP(1:list%NUMBER_IN_LIST)=list%LIST_SP(1:list%NUMBER_IN_LIST)
                      CALL MOVE_ALLOC(newListSP,list%LIST_SP)
                      list%SIZE=newSize
                    END IF
                    list%LIST_SP(list%NUMBER_IN_LIST+1:list%NUMBER_IN_LIST+appendedList%NUMBER_IN_LIST)= &
                      & appendedList%LIST_SP(1:appendedList%NUMBER_IN_LIST)
                    list%NUMBER_IN_LIST=list%NUMBER_IN_LIST+appendedList%NUMBER_IN_LIST
                  CASE(LIST_DP_TYPE)
                    IF(list%SIZE<list%NUMBER_IN_LIST+appendedList%NUMBER_IN_LIST) THEN
                      !Reallocate
                      newSize=MAX(2*list%NUMBER_IN_LIST,list%NUMBER_IN_LIST+appendedList%NUMBER_IN_LIST*2)
                      ALLOCATE(newListDP(newSize),stat=err)
                      IF(err/=0) CALL FLAG_ERROR("Could not allocate new list.",err,ERROR,*999)
                      newListDP(1:list%NUMBER_IN_LIST)=list%LIST_DP(1:list%NUMBER_IN_LIST)
                      CALL MOVE_ALLOC(newListDP,list%LIST_DP)
                      list%SIZE=newSize
                    END IF
                    list%LIST_DP(list%NUMBER_IN_LIST+1:list%NUMBER_IN_LIST+appendedList%NUMBER_IN_LIST)= &
                      & appendedList%LIST_DP(1:appendedList%NUMBER_IN_LIST)
                    list%NUMBER_IN_LIST=list%NUMBER_IN_LIST+appendedList%NUMBER_IN_LIST
                  CASE DEFAULT
                    CALL FLAG_ERROR("The list data type of "//TRIM(NUMBER_TO_VSTRING(list%DATA_TYPE,"*",err,error))// &
                      & " is invalid.",err,error,*999)
                  END SELECT
                CASE DEFAULT
                  CALL FLAG_ERROR("Dimensions > 1 not implemented for appended to a list",err,error,*999)
                END SELECT
              ELSE
                localError="Invalid data dimension. The list to append has data dimension of "// &
                  & TRIM(NUMBER_TO_VSTRING(appendedList%DATA_DIMENSION,"*",err,error))//" and the list data dimension is "// &
                  & TRIM(NUMBER_TO_VSTRING(list%DATA_DIMENSION,"*",err,error))//"."
                CALL FLAG_ERROR(localError,err,error,*999)
              ENDIF
            ELSE
              localError="The list data type of "//TRIM(NUMBER_TO_VSTRING(list%DATA_TYPE,"*",err,error))// &
                & " does not match the data type of the list to append"
              CALL FLAG_ERROR(localError,err,error,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("The list to append has not been finished",err,error,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("The list to append is not associated",err,error,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("The list has not been finished",err,error,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("List is not associated",err,error,*999)
    ENDIF

    CALL EXITS("List_AppendList")
    RETURN
999 IF(ALLOCATED(newListIntg)) DEALLOCATE(newListIntg)
    IF(ALLOCATED(newListSP)) DEALLOCATE(newListSP)
    IF(ALLOCATED(newListDP)) DEALLOCATE(newListDP)
    IF(ALLOCATED(newListCInt)) DEALLOCATE(newListCInt)
    CALL ERRORS("List_AppendList",err,error)
    CALL EXITS("List_AppendList")
    RETURN 1
  END SUBROUTINE List_AppendList

  !
  !================================================================================================================================
  !

  !>Clears all the items from a list
  SUBROUTINE List_ClearItems(list,err,error,*)
   !Argument Variables
    TYPE(LIST_TYPE), POINTER, INTENT(INOUT) :: list !<A pointer to the list
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string

    CALL ENTERS("List_ClearItems",err,error,*999)

    IF(ASSOCIATED(list)) THEN
      IF(list%LIST_FINISHED) THEN
        IF(list%mutable) THEN
          list%NUMBER_IN_LIST=0
        ELSE
          CALL FLAG_ERROR("The list is not mutable",err,error,*999)
        END IF
      ELSE
        CALL FLAG_ERROR("The list has not been finished",err,error,*999)
      END IF
    ELSE
      CALL FLAG_ERROR("List is not associated",err,error,*999)
    END IF

    CALL EXITS("List_ClearItems")
    RETURN
999 CALL ERRORS("List_ClearItems",err,error)
    CALL EXITS("List_ClearItems")
    RETURN 1
  END SUBROUTINE List_ClearItems
  !
  !================================================================================================================================
  !

  !>Initialises a list and all its components
  SUBROUTINE List_Initialise(list,err,error,*)

    !Argument Variables
    TYPE(LIST_TYPE), POINTER, INTENT(INOUT) :: list !<A pointer to the list to initialise. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError    

    CALL Enters("List_Initialise",err,error,*998)

    IF(ASSOCIATED(list)) THEN
      CALL FlagError("List is already associated.",err,error,*998)
    ELSE
      ALLOCATE(list,STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate list.",err,error,*999)
      list%LIST_FINISHED=.FALSE.
      list%mutable=.FALSE.
      list%NUMBER_IN_LIST=0
      list%DATA_DIMENSION=1
      list%INITIAL_SIZE=10
      list%size=0
      list%DATA_TYPE=LIST_INTG_TYPE
      list%KEY_DIMENSION=1
      list%SORT_ORDER=LIST_SORT_ASCENDING_TYPE
      list%SORT_METHOD=LIST_HEAP_SORT_METHOD
    ENDIF

    Call Exits("List_Initialise")
    RETURN
999 CALL List_Finalise(list,dummyErr,dummyError,*998)
998 CALL Errors("List_Initialise",err,error)
    Call Exits("List_Initialise")
    RETURN 1
  END SUBROUTINE List_Initialise

  !
  !================================================================================================================================
  !

  !>Sets/changes the initial size for a list
  SUBROUTINE List_InitialSizeSet(list,initialSize,err,error,*)

    !Argument Variables
    TYPE(LIST_TYPE), POINTER, INTENT(IN) :: list !<A pointer to the list
    INTEGER(INTG), INTENT(IN) :: initialSize !<The initial size of the list to set. Must be greater than zero.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    CALL Enters("List_InitialSizeSet",err,error,*999)

    IF(ASSOCIATED(list)) THEN
      IF(list%LIST_FINISHED) THEN
        CALL FlagError("List has been finished.",err,error,*999)
      ELSE
        IF(initialSize>0) THEN
          list%INITIAL_SIZE=initialSize
        ELSE
          localError="The initial size of "//TRIM(NumberToVString(initialSize,"*",err,error))// &
            & " is invalid. The size must be > 0."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("List is not associated",err,error,*999)
    ENDIF

    Call Exits("List_InitialSizeSet")
    RETURN
999 CALL Errors("List_InitialSizeSet",err,error)
    Call Exits("List_InitialSizeSet")
    RETURN 1
  END SUBROUTINE List_InitialSizeSet

  !
  !================================================================================================================================
  !

  !>Sets/changes the initial size for a list
  SUBROUTINE LIST_INITIAL_SIZE_SET(LIST,INITIAL_SIZE,err,error,*)

    !Argument Variables
    TYPE(LIST_TYPE), POINTER, INTENT(IN) :: LIST !<A pointer to the list
    INTEGER(INTG), INTENT(IN) :: INITIAL_SIZE !<The initial size of the list to set. Must be greater than zero.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL List_InitialSizeSet(LIST,INITIAL_SIZE,err,error,*999)
    
    RETURN
999 RETURN 1
  END SUBROUTINE LIST_INITIAL_SIZE_SET

  !
  !================================================================================================================================
  !

  !>Adds an item to the end of an integer list of data dimension 1. 
  SUBROUTINE LIST_ITEM_ADD_INTG1(LIST,ITEM,err,error,*)
   !Argument Variables
    TYPE(LIST_TYPE), POINTER, INTENT(INOUT) :: LIST !<A pointer to the list
    INTEGER(INTG), INTENT(IN) :: ITEM !<The item to add
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: NEW_SIZE
    INTEGER(INTG), ALLOCATABLE :: NEW_LIST(:)
    TYPE(VARYING_STRING) :: localError

    CALL Enters("LIST_ITEM_ADD_INTG1",err,error,*999)

    IF(ASSOCIATED(LIST)) THEN
      IF(list%LIST_FINISHED) THEN
        IF(list%DATA_TYPE==LIST_INTG_TYPE) THEN
          IF(list%DATA_DIMENSION==1) THEN
            IF(list%NUMBER_IN_LIST==list%SIZE) THEN
              !Reallocate
              NEW_SIZE=MAX(2*list%NUMBER_IN_LIST,1)
              ALLOCATE(NEW_LIST(NEW_SIZE),STAT=err)
              IF(err/=0) CALL FlagError("Could not allocate new list.",err,error,*999)
              NEW_LIST(1:list%NUMBER_IN_LIST)=list%LIST_INTG(1:list%NUMBER_IN_LIST)
              CALL MOVE_ALLOC(NEW_LIST,list%LIST_INTG)
              list%SIZE=NEW_SIZE
            ENDIF
            list%LIST_INTG(list%NUMBER_IN_LIST+1)=ITEM
            list%NUMBER_IN_LIST=list%NUMBER_IN_LIST+1
          ELSE
            localError="Invalid data dimension. The supplied data dimension is 1 and the list data dimension is "// &
              & TRIM(NumberToVString(list%DATA_DIMENSION,"*",err,error))//"."
            CALL FlagError(localError,err,error,*999)
          ENDIF
        ELSE
          localError="The list data type of "//TRIM(NumberToVString(list%DATA_TYPE,"*",err,error))// &
            & " does not match the integer type of the supplied list item"
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("The list has not been finished",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("List is not associated",err,error,*999)
    ENDIF
    
    Call Exits("LIST_ITEM_ADD_INTG1")
    RETURN
999 IF(ALLOCATED(NEW_LIST)) DEALLOCATE(NEW_LIST)
    CALL Errors("LIST_ITEM_ADD_INTG1",err,error)
    Call Exits("LIST_ITEM_ADD_INTG1")
    RETURN 1
  END SUBROUTINE LIST_ITEM_ADD_INTG1
  
  !
  !================================================================================================================================
  !

  !>Adds an item to the end of an integer list of data dimension > 1. 
  SUBROUTINE LIST_ITEM_ADD_INTG2(LIST,ITEM,err,error,*)
   !Argument Variables
    TYPE(LIST_TYPE), POINTER, INTENT(INOUT) :: LIST !<A pointer to the list
    INTEGER(INTG), INTENT(IN) :: ITEM(:) !<The item to add
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: NEW_SIZE
    INTEGER(INTG), ALLOCATABLE :: NEW_LIST(:,:)
    TYPE(VARYING_STRING) :: localError

    CALL Enters("LIST_ITEM_ADD_INTG2",err,error,*999)

    IF(ASSOCIATED(LIST)) THEN
      IF(list%LIST_FINISHED) THEN
        IF(list%DATA_TYPE==LIST_INTG_TYPE) THEN
          IF(list%DATA_DIMENSION==SIZE(ITEM,1)) THEN
            IF(list%NUMBER_IN_LIST==list%SIZE) THEN
              !Reallocate
              NEW_SIZE=MAX(2*list%NUMBER_IN_LIST,1)
              ALLOCATE(NEW_LIST(list%DATA_DIMENSION,NEW_SIZE),STAT=err)
              IF(err/=0) CALL FlagError("Could not allocate new list.",err,error,*999)
              NEW_LIST(:,1:list%NUMBER_IN_LIST)=list%LIST_INTG2(:,1:list%NUMBER_IN_LIST)
              CALL MOVE_ALLOC(NEW_LIST,list%LIST_INTG2)
              list%SIZE=NEW_SIZE
            ENDIF
            list%LIST_INTG2(:,list%NUMBER_IN_LIST+1)=ITEM
            list%NUMBER_IN_LIST=list%NUMBER_IN_LIST+1
          ELSE
            localError="Invalid data dimension. The supplied data dimension is "// &
              & TRIM(NumberToVString(SIZE(ITEM,1),"*",err,error))//" and the list data dimension is "// &
              & TRIM(NumberToVString(list%DATA_DIMENSION,"*",err,error))//"."
            CALL FlagError(localError,err,error,*999)
          ENDIF
        ELSE
          localError="The list data type of "//TRIM(NumberToVString(list%DATA_TYPE,"*",err,error))// &
            & " does not match the integer type of the supplied list item."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("The list has not been finished.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("List is not associated.",err,error,*999)
    ENDIF
    
    Call Exits("LIST_ITEM_ADD_INTG2")
    RETURN
999 IF(ALLOCATED(NEW_LIST)) DEALLOCATE(NEW_LIST)
    CALL Errors("LIST_ITEM_ADD_INTG2",err,error)
    Call Exits("LIST_ITEM_ADD_INTG2")
    RETURN 1
    
  END SUBROUTINE LIST_ITEM_ADD_INTG2
  
  !
  !================================================================================================================================
  !

  !>Adds an item to the end of a single precision real list of data dimension 1. 
  SUBROUTINE LIST_ITEM_ADD_SP1(LIST,ITEM,err,error,*)

    !Argument Variables
    TYPE(LIST_TYPE), POINTER, INTENT(INOUT) :: LIST !<A pointer to the list
    REAL(SP), INTENT(IN) :: ITEM !<The item to add
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: NEW_SIZE
    REAL(SP), ALLOCATABLE :: NEW_LIST(:)
    TYPE(VARYING_STRING) :: localError
    
    CALL Enters("LIST_ITEM_ADD_SP1",err,error,*999)

    IF(ASSOCIATED(LIST)) THEN
      IF(list%LIST_FINISHED) THEN
        IF(list%DATA_TYPE==LIST_SP_TYPE) THEN
          IF(list%DATA_DIMENSION==1) THEN
            IF(list%NUMBER_IN_LIST==list%SIZE) THEN
              !Reallocate
              NEW_SIZE=MAX(2*list%NUMBER_IN_LIST,1)
              ALLOCATE(NEW_LIST(NEW_SIZE),STAT=err)
              IF(err/=0) CALL FlagError("Could not allocate new list.",err,error,*999)
              NEW_LIST(1:list%NUMBER_IN_LIST)=list%LIST_SP(1:list%NUMBER_IN_LIST)
              CALL MOVE_ALLOC(NEW_LIST,list%LIST_SP)
              list%SIZE=NEW_SIZE
            ENDIF
            list%LIST_SP(list%NUMBER_IN_LIST+1)=ITEM
            list%NUMBER_IN_LIST=list%NUMBER_IN_LIST+1
          ELSE
            localError="Invalid data dimension. The supplied data dimension is 1 and the list data dimension is "// &
              & TRIM(NumberToVString(list%DATA_DIMENSION,"*",err,error))//"."
          ENDIF
        ELSE
          localError="The list data type of "//TRIM(NumberToVString(list%DATA_TYPE,"*",err,error))// &
            & " does not match the single precision type of the supplied list item."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("The list has not been finished.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("List is not associated.",err,error,*999)
    ENDIF
    Call Exits("LIST_ITEM_ADD_SP1")
    RETURN
999 IF(ALLOCATED(NEW_LIST)) DEALLOCATE(NEW_LIST)
    CALL Errors("LIST_ITEM_ADD_SP1",err,error)
    Call Exits("LIST_ITEM_ADD_SP1")
    RETURN 1
  END SUBROUTINE LIST_ITEM_ADD_SP1
  
  !
  !================================================================================================================================
  !

  !>Adds an item to the end of a single precision real list of data dimension > 1. 
  SUBROUTINE LIST_ITEM_ADD_SP2(LIST,ITEM,err,error,*)

    !Argument Variables
    TYPE(LIST_TYPE), POINTER, INTENT(INOUT) :: LIST !<A pointer to the list
    REAL(SP), INTENT(IN) :: ITEM(:) !<The item to add
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: NEW_SIZE
    REAL(SP), ALLOCATABLE :: NEW_LIST(:,:)
    TYPE(VARYING_STRING) :: localError
    
    CALL Enters("LIST_ITEM_ADD_SP2",err,error,*999)

    IF(ASSOCIATED(LIST)) THEN
      IF(list%LIST_FINISHED) THEN
        IF(list%DATA_TYPE==LIST_SP_TYPE) THEN
          IF(list%DATA_DIMENSION==SIZE(ITEM,1)) THEN
            IF(list%NUMBER_IN_LIST==list%SIZE) THEN
              !Reallocate
              NEW_SIZE=MAX(2*list%NUMBER_IN_LIST,1)
              ALLOCATE(NEW_LIST(list%DATA_DIMENSION,NEW_SIZE),STAT=err)
              IF(err/=0) CALL FlagError("Could not allocate new list.",err,error,*999)
              NEW_LIST(:,1:list%NUMBER_IN_LIST)=list%LIST_SP2(:,1:list%NUMBER_IN_LIST)
              CALL MOVE_ALLOC(NEW_LIST,list%LIST_SP2)
              list%SIZE=NEW_SIZE
            ENDIF
            list%LIST_SP2(:,list%NUMBER_IN_LIST+1)=ITEM
            list%NUMBER_IN_LIST=list%NUMBER_IN_LIST+1
          ELSE
            localError="Invalid data dimension. The supplied data dimension is "// &
              & TRIM(NumberToVString(SIZE(ITEM,1),"*",err,error))//" and the list data dimension is "// &
              & TRIM(NumberToVString(list%DATA_DIMENSION,"*",err,error))//"."
            CALL FlagError(localError,err,error,*999)
          ENDIF
        ELSE
          localError="The list data type of "//TRIM(NumberToVString(list%DATA_TYPE,"*",err,error))// &
            & " does not match the single precision type of the supplied list item."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("The list has not been finished.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("List is not associated.",err,error,*999)
    ENDIF
    Call Exits("LIST_ITEM_ADD_SP2")
    RETURN
999 IF(ALLOCATED(NEW_LIST)) DEALLOCATE(NEW_LIST)
    CALL Errors("LIST_ITEM_ADD_SP2",err,error)
    Call Exits("LIST_ITEM_ADD_SP2")
    RETURN 1
  END SUBROUTINE LIST_ITEM_ADD_SP2
  
  !
  !================================================================================================================================
  !

  !>Adds an item to the end of a double precision real list of data dimension 1.
  SUBROUTINE LIST_ITEM_ADD_DP1(LIST,ITEM,err,error,*)

    !Argument Variables
    TYPE(LIST_TYPE), POINTER, INTENT(INOUT) :: LIST !<A pointer to the list
    REAL(DP), INTENT(IN) :: ITEM !<The item to add
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: NEW_SIZE
    REAL(DP), ALLOCATABLE :: NEW_LIST(:)
    TYPE(VARYING_STRING) :: localError

    CALL Enters("LIST_ITEM_ADD_DP1",err,error,*999)

    IF(ASSOCIATED(LIST)) THEN
      IF(list%LIST_FINISHED) THEN
        IF(list%DATA_TYPE==LIST_DP_TYPE) THEN
          IF(list%DATA_DIMENSION==1) THEN
            IF(list%NUMBER_IN_LIST==list%SIZE) THEN
              !Reallocate
              NEW_SIZE=MAX(2*list%NUMBER_IN_LIST,1)
              ALLOCATE(NEW_LIST(NEW_SIZE),STAT=err)
              IF(err/=0) CALL FlagError("Could not allocate new list.",err,error,*999)
              NEW_LIST(1:list%NUMBER_IN_LIST)=list%LIST_DP(1:list%NUMBER_IN_LIST)
              CALL MOVE_ALLOC(NEW_LIST,list%LIST_DP)
              list%SIZE=NEW_SIZE
            ENDIF
            list%LIST_DP(list%NUMBER_IN_LIST+1)=ITEM
            list%NUMBER_IN_LIST=list%NUMBER_IN_LIST+1
          ELSE
            localError="Invalid data dimension. The supplied data dimension is 1 and the list data dimension is "// &
              & TRIM(NumberToVString(list%DATA_DIMENSION,"*",err,error))//"."
            CALL FlagError(localError,err,error,*999)
          ENDIF
        ELSE
          localError="The list data type of "//TRIM(NumberToVString(list%DATA_TYPE,"*",err,error))// &
            & " does not match the double precision type of the supplied list item."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("The list has not been finished.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("List is not associated.",err,error,*999)
    ENDIF
    Call Exits("LIST_ITEM_ADD_DP1")
    RETURN
999 IF(ALLOCATED(NEW_LIST)) DEALLOCATE(NEW_LIST)
    CALL Errors("LIST_ITEM_ADD_DP1",err,error)
    Call Exits("LIST_ITEM_ADD_DP1")
    RETURN 1
  END SUBROUTINE LIST_ITEM_ADD_DP1
  
  !
  !================================================================================================================================
  !

  !>Adds an item to the end of a double precision real list of data dimension > 1.
  SUBROUTINE LIST_ITEM_ADD_DP2(LIST,ITEM,err,error,*)

    !Argument Variables
    TYPE(LIST_TYPE), POINTER, INTENT(INOUT) :: LIST !<A pointer to the list
    REAL(DP), INTENT(IN) :: ITEM(:) !<The item to add
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: NEW_SIZE
    REAL(DP), ALLOCATABLE :: NEW_LIST(:,:)
    TYPE(VARYING_STRING) :: localError

    CALL Enters("LIST_ITEM_ADD_DP2",err,error,*999)

    IF(ASSOCIATED(LIST)) THEN
      IF(list%LIST_FINISHED) THEN
        IF(list%DATA_TYPE==LIST_DP_TYPE) THEN
          IF(list%DATA_DIMENSION==SIZE(ITEM,1)) THEN
            IF(list%NUMBER_IN_LIST==list%SIZE) THEN
              !Reallocate
              NEW_SIZE=MAX(2*list%NUMBER_IN_LIST,1)
              ALLOCATE(NEW_LIST(list%DATA_DIMENSION,NEW_SIZE),STAT=err)
              IF(err/=0) CALL FlagError("Could not allocate new list.",err,error,*999)
              NEW_LIST(:,1:list%NUMBER_IN_LIST)=list%LIST_DP2(:,1:list%NUMBER_IN_LIST)
              CALL MOVE_ALLOC(NEW_LIST,list%LIST_DP2)
              list%SIZE=NEW_SIZE
            ENDIF
            list%LIST_DP2(:,list%NUMBER_IN_LIST+1)=ITEM
            list%NUMBER_IN_LIST=list%NUMBER_IN_LIST+1
          ELSE
            localError="Invalid data dimension. The supplied data dimension is "// &
              & TRIM(NumberToVString(SIZE(ITEM,1),"*",err,error))//" and the list data dimension is "// &
              & TRIM(NumberToVString(list%DATA_DIMENSION,"*",err,error))//"."
            CALL FlagError(localError,err,error,*999)
          ENDIF
        ELSE
          localError="The list data type of "//TRIM(NumberToVString(list%DATA_TYPE,"*",err,error))// &
            & " does not match the double precision type of the supplied list item."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("The list has not been finished.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("List is not associated.",err,error,*999)
    ENDIF
    Call Exits("LIST_ITEM_ADD_DP2")
    RETURN
999 IF(ALLOCATED(NEW_LIST)) DEALLOCATE(NEW_LIST)
    CALL Errors("LIST_ITEM_ADD_DP2",err,error)
    Call Exits("LIST_ITEM_ADD_DP2")
    RETURN 1
  END SUBROUTINE LIST_ITEM_ADD_DP2
  
  !
  !================================================================================================================================
  !

  !>Sets an item in an integer list of data dimension 1. 
  SUBROUTINE LIST_ITEM_SET_INTG1(LIST,LIST_ITEM,ITEM,err,error,*)
   !Argument Variables
    TYPE(LIST_TYPE), POINTER, INTENT(IN) :: LIST !<A pointer to the list
    INTEGER(INTG), INTENT(IN) :: LIST_ITEM !<The index of the item to set
    INTEGER(INTG), INTENT(IN) :: ITEM !<The item to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    CALL Enters("LIST_ITEM_SET_INTG1",err,error,*999)

    IF(ASSOCIATED(LIST)) THEN
      IF(list%LIST_FINISHED) THEN
        IF(list%DATA_TYPE==LIST_INTG_TYPE) THEN
          IF(list%DATA_DIMENSION==1) THEN
            IF(LIST_ITEM>0.AND.LIST_ITEM<=list%NUMBER_IN_LIST) THEN
              IF(list%MUTABLE) THEN
                list%LIST_INTG(LIST_ITEM)=ITEM
              ELSE
                CALL FlagError("Cannot modify an immutable list.",err,error,*999)
              ENDIF
            ELSE
              localError="Invalid list index. The supplied index is "// &
                & TRIM(NumberToVString(LIST_ITEM,"*",err,error))//" and that list entry count is"// &
                & TRIM(NumberToVString(list%NUMBER_IN_LIST,"*",err,error))//"."
              CALL FlagError(localError,err,error,*999)
            ENDIF
          ELSE
            localError="Invalid data dimension. The supplied data dimension is 1 and the list data dimension is "// &
              & TRIM(NumberToVString(list%DATA_DIMENSION,"*",err,error))//"."
            CALL FlagError(localError,err,error,*999)
          ENDIF
        ELSE
          localError="The list data type of "//TRIM(NumberToVString(list%DATA_TYPE,"*",err,error))// &
            & " does not match the integer type of the supplied list item"
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("The list has not been finished",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("List is not associated",err,error,*999)
    ENDIF
    
    Call Exits("LIST_ITEM_SET_INTG1")
    RETURN
999 CALL Errors("LIST_ITEM_SET_INTG1",err,error)
    Call Exits("LIST_ITEM_SET_INTG1")
    RETURN 1
  END SUBROUTINE LIST_ITEM_SET_INTG1
  
  !
  !================================================================================================================================
  !

  !>Set an item in an integer list of data dimension > 1. 
  SUBROUTINE LIST_ITEM_SET_INTG2(LIST,LIST_ITEM,ITEM,err,error,*)
   !Argument Variables
    TYPE(LIST_TYPE), POINTER, INTENT(IN) :: LIST !<A pointer to the list
    INTEGER(INTG), INTENT(IN) :: LIST_ITEM !<The index of the item to set
    INTEGER(INTG), INTENT(IN) :: ITEM(:) !<The item to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    CALL Enters("LIST_ITEM_SET_INTG2",err,error,*999)

    IF(ASSOCIATED(LIST)) THEN
      IF(list%LIST_FINISHED) THEN
        IF(list%DATA_TYPE==LIST_INTG_TYPE) THEN
          IF(list%DATA_DIMENSION==SIZE(ITEM,1)) THEN
            IF(LIST_ITEM>0.AND.LIST_ITEM<=list%NUMBER_IN_LIST) THEN
              IF(list%MUTABLE) THEN
                list%LIST_INTG2(:,LIST_ITEM)=ITEM
              ELSE
                CALL FlagError("Cannot modify an immutable list.",err,error,*999)
              ENDIF
            ELSE
              localError="Invalid list index. The supplied index is "//&
                & TRIM(NumberToVString(LIST_ITEM,"*",err,error))//" and that list entry count is"//&
                & TRIM(NumberToVString(list%NUMBER_IN_LIST,"*",err,error))//"."
              CALL FlagError(localError,err,error,*999)
            ENDIF
          ELSE
            localError="Invalid data dimension. The supplied data dimension is "// &
              & TRIM(NumberToVString(SIZE(ITEM,1),"*",err,error))//" and the list data dimension is "// &
              & TRIM(NumberToVString(list%DATA_DIMENSION,"*",err,error))//"."
            CALL FlagError(localError,err,error,*999)
          ENDIF
        ELSE
          localError="The list data type of "//TRIM(NumberToVString(list%DATA_TYPE,"*",err,error))// &
            & " does not match the integer type of the supplied list item."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("The list has not been finished.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("List is not associated.",err,error,*999)
    ENDIF
    
    Call Exits("LIST_ITEM_SET_INTG2")
    RETURN
999 CALL Errors("LIST_ITEM_SET_INTG2",err,error)
    Call Exits("LIST_ITEM_SET_INTG2")
    RETURN 1
    
  END SUBROUTINE LIST_ITEM_SET_INTG2
  
  !
  !================================================================================================================================
  !

  !>Sets an item in a single precision real list of data dimension 1. 
  SUBROUTINE LIST_ITEM_SET_SP1(LIST,LIST_ITEM,ITEM,err,error,*)

    !Argument Variables
    TYPE(LIST_TYPE), POINTER, INTENT(IN) :: LIST !<A pointer to the list
    INTEGER(INTG), INTENT(IN) :: LIST_ITEM !<The index of the item to set
    REAL(SP), INTENT(IN) :: ITEM !<The item to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    CALL Enters("LIST_ITEM_SET_SP1",err,error,*999)

    IF(ASSOCIATED(LIST)) THEN
      IF(list%LIST_FINISHED) THEN
        IF(list%DATA_TYPE==LIST_SP_TYPE) THEN
          IF(list%DATA_DIMENSION==1) THEN
            IF(LIST_ITEM>0.AND.LIST_ITEM<=list%NUMBER_IN_LIST) THEN
              IF(list%MUTABLE) THEN
                list%LIST_SP(LIST_ITEM)=ITEM
              ELSE
                CALL FlagError("Cannot modify an immutable list.",err,error,*999)
              ENDIF
            ELSE
              localError="Invalid list index. The supplied index is "//&
                & TRIM(NumberToVString(LIST_ITEM,"*",err,error))//" and that list entry count is"//&
                & TRIM(NumberToVString(list%NUMBER_IN_LIST,"*",err,error))//"."
              CALL FlagError(localError,err,error,*999)
            ENDIF
          ELSE
            localError="Invalid data dimension. The supplied data dimension is 1 and the list data dimension is "// &
              & TRIM(NumberToVString(list%DATA_DIMENSION,"*",err,error))//"."
          ENDIF
        ELSE
          localError="The list data type of "//TRIM(NumberToVString(list%DATA_TYPE,"*",err,error))// &
            & " does not match the single precision type of the supplied list item."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("The list has not been finished.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("List is not associated.",err,error,*999)
    ENDIF
    Call Exits("LIST_ITEM_SET_SP1")
    RETURN
999 CALL Errors("LIST_ITEM_SET_SP1",err,error)
    Call Exits("LIST_ITEM_SET_SP1")
    RETURN 1
  END SUBROUTINE LIST_ITEM_SET_SP1
  
  !
  !================================================================================================================================
  !

  !>Sets an item in a single precision real list of data dimension > 1. 
  SUBROUTINE LIST_ITEM_SET_SP2(LIST,LIST_ITEM,ITEM,err,error,*)

    !Argument Variables
    TYPE(LIST_TYPE), POINTER, INTENT(IN) :: LIST !<A pointer to the list
    INTEGER(INTG), INTENT(IN) :: LIST_ITEM !<The index of the item to set
    REAL(SP), INTENT(IN) :: ITEM(:) !<The item to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    CALL Enters("LIST_ITEM_SET_SP2",err,error,*999)

    IF(ASSOCIATED(LIST)) THEN
      IF(list%LIST_FINISHED) THEN
        IF(list%DATA_TYPE==LIST_SP_TYPE) THEN
          IF(list%DATA_DIMENSION==SIZE(ITEM,1)) THEN
            IF(LIST_ITEM>0.AND.LIST_ITEM<=list%NUMBER_IN_LIST) THEN
              IF(list%MUTABLE) THEN
                list%LIST_SP2(:,LIST_ITEM)=ITEM
              ELSE
                CALL FlagError("Cannot modify an immutable list.",err,error,*999)
              ENDIF
            ELSE
              localError="Invalid list index. The supplied index is "//&
                & TRIM(NumberToVString(LIST_ITEM,"*",err,error))//" and that list entry count is"//&
                & TRIM(NumberToVString(list%NUMBER_IN_LIST,"*",err,error))//"."
              CALL FlagError(localError,err,error,*999)
            ENDIF
          ELSE
            localError="Invalid data dimension. The supplied data dimension is "// &
              & TRIM(NumberToVString(SIZE(ITEM,1),"*",err,error))//" and the list data dimension is "// &
              & TRIM(NumberToVString(list%DATA_DIMENSION,"*",err,error))//"."
            CALL FlagError(localError,err,error,*999)
          ENDIF
        ELSE
          localError="The list data type of "//TRIM(NumberToVString(list%DATA_TYPE,"*",err,error))// &
            & " does not match the single precision type of the supplied list item."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("The list has not been finished.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("List is not associated.",err,error,*999)
    ENDIF
    Call Exits("LIST_ITEM_SET_SP2")
    RETURN
999 CALL Errors("LIST_ITEM_SET_SP2",err,error)
    Call Exits("LIST_ITEM_SET_SP2")
    RETURN 1
  END SUBROUTINE LIST_ITEM_SET_SP2
  
  !
  !================================================================================================================================
  !

  !>Sets an item in a double precision real list of data dimension 1.
  SUBROUTINE LIST_ITEM_SET_DP1(LIST,LIST_ITEM,ITEM,err,error,*)

    !Argument Variables
    TYPE(LIST_TYPE), POINTER, INTENT(IN) :: LIST !<A pointer to the list
    INTEGER(INTG), INTENT(IN) :: LIST_ITEM !<The index of the item to set. 
    REAL(DP), INTENT(IN) :: ITEM !<The item to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    CALL Enters("LIST_ITEM_SET_DP1",err,error,*999)

    IF(ASSOCIATED(LIST)) THEN
      IF(list%LIST_FINISHED) THEN
        IF(list%DATA_TYPE==LIST_DP_TYPE) THEN
          IF(list%DATA_DIMENSION==1) THEN
            IF(LIST_ITEM>0.AND.LIST_ITEM<=list%NUMBER_IN_LIST) THEN
              IF(list%MUTABLE) THEN
                list%LIST_DP(LIST_ITEM)=ITEM
              ELSE
                CALL FlagError("Cannot modify an immutable list.",err,error,*999)
              ENDIF
            ELSE
              localError="Invalid list index. The supplied index is "//&
                & TRIM(NumberToVString(LIST_ITEM,"*",err,error))//" and that list entry count is"//&
                & TRIM(NumberToVString(list%NUMBER_IN_LIST,"*",err,error))//"."
              CALL FlagError(localError,err,error,*999)
            ENDIF
          ELSE
            localError="Invalid data dimension. The supplied data dimension is 1 and the list data dimension is "// &
              & TRIM(NumberToVString(list%DATA_DIMENSION,"*",err,error))//"."
            CALL FlagError(localError,err,error,*999)
          ENDIF
        ELSE
          localError="The list data type of "//TRIM(NumberToVString(list%DATA_TYPE,"*",err,error))// &
            & " does not match the double precision type of the supplied list item."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("The list has not been finished.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("List is not associated.",err,error,*999)
    ENDIF
    Call Exits("LIST_ITEM_SET_DP1")
    RETURN
999 CALL Errors("LIST_ITEM_SET_DP1",err,error)
    Call Exits("LIST_ITEM_SET_DP1")
    RETURN 1
  END SUBROUTINE LIST_ITEM_SET_DP1
  
  !
  !================================================================================================================================
  !

  !>Sets an item in a double precision real list of data dimension > 1.
  SUBROUTINE LIST_ITEM_SET_DP2(LIST,LIST_ITEM,ITEM,err,error,*)

    !Argument Variables
    TYPE(LIST_TYPE), POINTER, INTENT(IN) :: LIST !<A pointer to the list
    INTEGER(INTG), INTENT(IN) :: LIST_ITEM !<The index of the item to set.
    REAL(DP), INTENT(IN) :: ITEM(:) !<The item to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    CALL Enters("LIST_ITEM_SET_DP2",err,error,*999)

    IF(ASSOCIATED(LIST)) THEN
      IF(list%LIST_FINISHED) THEN
        IF(list%DATA_TYPE==LIST_DP_TYPE) THEN
          IF(list%DATA_DIMENSION==SIZE(ITEM,1)) THEN
            IF(LIST_ITEM>0.AND.LIST_ITEM<=list%NUMBER_IN_LIST) THEN
              IF(list%MUTABLE) THEN
                list%LIST_DP2(:,LIST_ITEM)=ITEM
              ELSE
                CALL FlagError("Cannot modify an immutable list.",err,error,*999)
              ENDIF
            ELSE
              localError="Invalid list index. The supplied index is "//&
                & TRIM(NumberToVString(LIST_ITEM,"*",err,error))//" and that list entry count is"//&
                & TRIM(NumberToVString(list%NUMBER_IN_LIST,"*",err,error))//"."
              CALL FlagError(localError,err,error,*999)
            ENDIF
          ELSE
            localError="Invalid data dimension. The supplied data dimension is "// &
              & TRIM(NumberToVString(SIZE(ITEM,1),"*",err,error))//" and the list data dimension is "// &
              & TRIM(NumberToVString(list%DATA_DIMENSION,"*",err,error))//"."
            CALL FlagError(localError,err,error,*999)
          ENDIF
        ELSE
          localError="The list data type of "//TRIM(NumberToVString(list%DATA_TYPE,"*",err,error))// &
            & " does not match the double precision type of the supplied list item."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("The list has not been finished.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("List is not associated.",err,error,*999)
    ENDIF
    Call Exits("LIST_ITEM_SET_DP2")
    RETURN
999 CALL Errors("LIST_ITEM_SET_DP2",err,error)
    Call Exits("LIST_ITEM_SET_DP2")
    RETURN 1
  END SUBROUTINE LIST_ITEM_SET_DP2
  
  !
  !================================================================================================================================
  !

  !>Returns the ITEM in a list at position LIST_ITEM in the given integer LIST. 
  SUBROUTINE LIST_ITEM_GET_INTG1(LIST,LIST_ITEM,ITEM,err,error,*)

    !Argument Variables
    TYPE(LIST_TYPE), POINTER, INTENT(IN) :: LIST !<The pointer to the list
    INTEGER(INTG), INTENT(IN) :: LIST_ITEM !<The position of the item to get
    INTEGER(INTG), INTENT(OUT) :: ITEM !<On exit, the item at the specified position
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    CALL Enters("LIST_ITEM_GET_INTG1",err,error,*999)

    IF(ASSOCIATED(LIST)) THEN
      IF(list%LIST_FINISHED) THEN
        IF(list%DATA_TYPE==LIST_INTG_TYPE) THEN
          IF(LIST_ITEM>0.AND.LIST_ITEM<=list%NUMBER_IN_LIST) THEN
            IF(list%DATA_DIMENSION==1) THEN
              ITEM=list%LIST_INTG(LIST_ITEM)
            ELSE
              localError="Invalid item dimension. The specified item has dimension 1 and the list is of dimension "// &
                & TRIM(NumberToVString(list%DATA_DIMENSION,"*",err,error))
              CALL FlagError(localError,err,error,*999)
            ENDIF
          ELSE
            localError="The specified list item position of "//TRIM(NumberToVString(LIST_ITEM,"*",err,error))// &
              & " is invalid. The list item position must be > 0 and <= "// &
              & TRIM(NumberToVString(list%NUMBER_IN_LIST,"*",err,error))//"."
            CALL FlagError(localError,err,error,*999)
          ENDIF
        ELSE
          localError="The list data type of "//TRIM(NumberToVString(list%DATA_TYPE,"*",err,error))// &
            & " does not match the integer type of the supplied list item."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("List has not been finished.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("List is not associated.",err,error,*999)
    ENDIF

    Call Exits("LIST_ITEM_GET_INTG1")
    RETURN
999 CALL Errors("LIST_ITEM_GET_INTG1",err,error)
    Call Exits("LIST_ITEM_GET_INTG1")
    RETURN 1
  END SUBROUTINE LIST_ITEM_GET_INTG1
  
  !
  !================================================================================================================================
  !

  !>Returns the ITEM in a list at position LIST_ITEM in the given integer LIST. 
  SUBROUTINE LIST_ITEM_GET_INTG2(LIST,LIST_ITEM,ITEM,err,error,*)

    !Argument Variables
    TYPE(LIST_TYPE), POINTER, INTENT(IN) :: LIST !<The pointer to the list
    INTEGER(INTG), INTENT(IN) :: LIST_ITEM !<The position of the item to get
    INTEGER(INTG), INTENT(OUT) :: ITEM(:) !<On exit, the item at the specified position
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    CALL Enters("LIST_ITEM_GET_INTG2",err,error,*999)

    IF(ASSOCIATED(LIST)) THEN
      IF(list%LIST_FINISHED) THEN
        IF(list%DATA_TYPE==LIST_INTG_TYPE) THEN
          IF(LIST_ITEM>0.AND.LIST_ITEM<=list%NUMBER_IN_LIST) THEN
            IF(list%DATA_DIMENSION==SIZE(ITEM,1)) THEN
              ITEM=list%LIST_INTG2(:,LIST_ITEM)
            ELSE
              localError="Invalid item dimension. The specified item has dimension "// &
                & TRIM(NumberToVString(SIZE(ITEM,1),"*",err,error))//" and the list is of dimension "// &
                & TRIM(NumberToVString(list%DATA_DIMENSION,"*",err,error))
              CALL FlagError(localError,err,error,*999)
            ENDIF
          ELSE
            localError="The specified list item position of "//TRIM(NumberToVString(LIST_ITEM,"*",err,error))// &
              & " is invalid. The list item position must be > 0 and <= "// &
              & TRIM(NumberToVString(list%NUMBER_IN_LIST,"*",err,error))//"."
            CALL FlagError(localError,err,error,*999)
          ENDIF
        ELSE
          localError="The list data type of "//TRIM(NumberToVString(list%DATA_TYPE,"*",err,error))// &
            & " does not match the integer type of the supplied list item."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("List has not been finished.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("List is not associated.",err,error,*999)
    ENDIF

    Call Exits("LIST_ITEM_GET_INTG2")
    RETURN
999 CALL Errors("LIST_ITEM_GET_INTG2",err,error)
    Call Exits("LIST_ITEM_GET_INTG2")
    RETURN 1
  END SUBROUTINE LIST_ITEM_GET_INTG2
  
  !
  !================================================================================================================================
  !

  !>Returns the ITEM in a list at position LIST_ITEM in the given single precision LIST. 
  SUBROUTINE LIST_ITEM_GET_SP1(LIST,LIST_ITEM,ITEM,err,error,*)

    !Argument Variables
    TYPE(LIST_TYPE), POINTER, INTENT(IN) :: LIST !<The pointer to the list
    INTEGER(INTG), INTENT(IN) :: LIST_ITEM !<The position of the item to get
    REAL(SP), INTENT(OUT) :: ITEM !<On exit, the item at the specified position.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    CALL Enters("LIST_ITEM_GET_SP1",err,error,*999)

    IF(ASSOCIATED(LIST)) THEN
      IF(list%LIST_FINISHED) THEN
        IF(list%DATA_TYPE==LIST_SP_TYPE) THEN
          IF(LIST_ITEM>0.AND.LIST_ITEM<=list%NUMBER_IN_LIST) THEN
            IF(list%DATA_DIMENSION==1) THEN
              ITEM=list%LIST_SP(LIST_ITEM)
            ELSE
              localError="Invalid item dimension. The specified item has dimension 1 and the list is of dimension "// &
                & TRIM(NumberToVString(list%DATA_DIMENSION,"*",err,error))
              CALL FlagError(localError,err,error,*999)
            ENDIF
          ELSE
            localError="The specified list item position of "//TRIM(NumberToVString(LIST_ITEM,"*",err,error))// &
              & " is invalid. The list item position must be > 0 and <= "// &
              & TRIM(NumberToVString(list%NUMBER_IN_LIST,"*",err,error))//"."
            CALL FlagError(localError,err,error,*999)
          ENDIF
        ELSE
          localError="The list data type of "//TRIM(NumberToVString(list%DATA_TYPE,"*",err,error))// &
            & " does not match the single precision type of the supplied list item."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("List has not been finished.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("List is not associated.",err,error,*999)
    ENDIF

    Call Exits("LIST_ITEM_GET_SP1")
    RETURN
999 CALL Errors("LIST_ITEM_GET_SP1",err,error)
    Call Exits("LIST_ITEM_GET_SP1")
    RETURN 1
  END SUBROUTINE LIST_ITEM_GET_SP1
  
  !
  !================================================================================================================================
  !

  !>Returns the ITEM in a list at position LIST_ITEM in the given single precision LIST. 
  SUBROUTINE LIST_ITEM_GET_SP2(LIST,LIST_ITEM,ITEM,err,error,*)

    !Argument Variables
    TYPE(LIST_TYPE), POINTER, INTENT(IN) :: LIST !<The pointer to the list
    INTEGER(INTG), INTENT(IN) :: LIST_ITEM !<The position of the item to get
    REAL(SP), INTENT(OUT) :: ITEM(:) !<On exit, the item at the specified position
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    CALL Enters("LIST_ITEM_GET_SP2",err,error,*999)

    IF(ASSOCIATED(LIST)) THEN
      IF(list%LIST_FINISHED) THEN
        IF(list%DATA_TYPE==LIST_SP_TYPE) THEN
          IF(LIST_ITEM>0.AND.LIST_ITEM<=list%NUMBER_IN_LIST) THEN
            IF(list%DATA_DIMENSION==SIZE(ITEM,1)) THEN
              ITEM=list%LIST_SP2(:,LIST_ITEM)
            ELSE
              localError="Invalid item dimension. The specified item has dimension "// &
                & TRIM(NumberToVString(SIZE(ITEM,1),"*",err,error))//" and the list is of dimension "// &
                & TRIM(NumberToVString(list%DATA_DIMENSION,"*",err,error))
              CALL FlagError(localError,err,error,*999)
            ENDIF
          ELSE
            localError="The specified list item position of "//TRIM(NumberToVString(LIST_ITEM,"*",err,error))// &
              & " is invalid. The list item position must be > 0 and <= "// &
              & TRIM(NumberToVString(list%NUMBER_IN_LIST,"*",err,error))//"."
            CALL FlagError(localError,err,error,*999)
          ENDIF
        ELSE
          localError="The list data type of "//TRIM(NumberToVString(list%DATA_TYPE,"*",err,error))// &
            & " does not match the single precision type of the supplied list item."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("List has not been finished.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("List is not associated.",err,error,*999)
    ENDIF

    Call Exits("LIST_ITEM_GET_SP2")
    RETURN
999 CALL Errors("LIST_ITEM_GET_SP2",err,error)
    Call Exits("LIST_ITEM_GET_SP2")
    RETURN 1
  END SUBROUTINE LIST_ITEM_GET_SP2
  
  !
  !================================================================================================================================
  !

  !>Returns the ITEM in a list at position LIST_ITEM in the given double precision LIST. 
  SUBROUTINE LIST_ITEM_GET_DP1(LIST,LIST_ITEM,ITEM,err,error,*)

    !Argument Variables
    TYPE(LIST_TYPE), POINTER, INTENT(IN) :: LIST !<The pointer to the list
    INTEGER(INTG), INTENT(IN) :: LIST_ITEM !<The position of the item to get
    REAL(DP), INTENT(OUT) :: ITEM !<On exit, the item at the specified position.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    CALL Enters("LIST_ITEM_GET_DP1",err,error,*999)

    IF(ASSOCIATED(LIST)) THEN
      IF(list%LIST_FINISHED) THEN
        IF(list%DATA_TYPE==LIST_DP_TYPE) THEN
          IF(LIST_ITEM>0.AND.LIST_ITEM<=list%NUMBER_IN_LIST) THEN
            IF(list%DATA_DIMENSION==1) THEN
              ITEM=list%LIST_DP(LIST_ITEM)
            ELSE
              localError="Invalid item dimension. The specified item has dimension 1 and the list is of dimension "// &
                & TRIM(NumberToVString(list%DATA_DIMENSION,"*",err,error))
              CALL FlagError(localError,err,error,*999)
            ENDIF
          ELSE
            localError="The specified list item position of "//TRIM(NumberToVString(LIST_ITEM,"*",err,error))// &
              & " is invalid. The list item position must be > 0 and <= "// &
              & TRIM(NumberToVString(list%NUMBER_IN_LIST,"*",err,error))//"."
            CALL FlagError(localError,err,error,*999)
          ENDIF
        ELSE
          localError="The list data type of "//TRIM(NumberToVString(list%DATA_TYPE,"*",err,error))// &
            & " does not match the double precision type of the supplied list item."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("List has not been finished.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("List is not associated.",err,error,*999)
    ENDIF

    Call Exits("LIST_ITEM_GET_DP1")
    RETURN
999 CALL Errors("LIST_ITEM_GET_DP1",err,error)
    Call Exits("LIST_ITEM_GET_DP1")
    RETURN 1
  END SUBROUTINE LIST_ITEM_GET_DP1
  
  !
  !================================================================================================================================
  !

  !>Returns the ITEM in a list at position LIST_ITEM in the given double precision LIST. 
  SUBROUTINE LIST_ITEM_GET_DP2(LIST,LIST_ITEM,ITEM,err,error,*)

    !Argument Variables
    TYPE(LIST_TYPE), POINTER, INTENT(IN) :: LIST !<The pointer to the list
    INTEGER(INTG), INTENT(IN) :: LIST_ITEM !<The position of the item to get
    REAL(DP), INTENT(OUT) :: ITEM(:) !<On exit, the item at the specified position
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    CALL Enters("LIST_ITEM_GET_DP2",err,error,*999)

    IF(ASSOCIATED(LIST)) THEN
      IF(list%LIST_FINISHED) THEN
        IF(list%DATA_TYPE==LIST_DP_TYPE) THEN
          IF(LIST_ITEM>0.AND.LIST_ITEM<=list%NUMBER_IN_LIST) THEN
            IF(list%DATA_DIMENSION==SIZE(ITEM,1)) THEN
              ITEM=list%LIST_DP2(:,LIST_ITEM)
            ELSE
              localError="Invalid item dimension. The specified item has dimension "// &
                & TRIM(NumberToVString(SIZE(ITEM,1),"*",err,error))//" and the list is of dimension "// &
                & TRIM(NumberToVString(list%DATA_DIMENSION,"*",err,error))
              CALL FlagError(localError,err,error,*999)
            ENDIF
          ELSE
            localError="The specified list item position of "//TRIM(NumberToVString(LIST_ITEM,"*",err,error))// &
              & " is invalid. The list item position must be > 0 and <= "// &
              & TRIM(NumberToVString(list%NUMBER_IN_LIST,"*",err,error))//"."
            CALL FlagError(localError,err,error,*999)
          ENDIF
        ELSE
          localError="The list data type of "//TRIM(NumberToVString(list%DATA_TYPE,"*",err,error))// &
            & " does not match the double precision type of the supplied list item."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("List has not been finished.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("List is not associated.",err,error,*999)
    ENDIF

    Call Exits("LIST_ITEM_GET_DP2")
    RETURN
999 CALL Errors("LIST_ITEM_GET_DP2",err,error)
    Call Exits("LIST_ITEM_GET_DP2")
    RETURN 1
  END SUBROUTINE LIST_ITEM_GET_DP2
  
  !
  !================================================================================================================================
  !

  !>Determines if ITEM is in the given integer LIST. If it is LIST_ITEM is the index in the list. If not LIST_ITEM is 0.
  SUBROUTINE LIST_ITEM_IN_LIST_INTG1(LIST,ITEM,LIST_ITEM,err,error,*)

    !Argument Variables
    TYPE(LIST_TYPE), POINTER, INTENT(IN) :: LIST !<The pointer to the list
    INTEGER(INTG), INTENT(IN) :: ITEM !<The item to find.
    INTEGER(INTG), INTENT(OUT) :: LIST_ITEM !<On exit, the position of the item in the list. If the item does not exist then the value of 0 is returned.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    CALL Enters("LIST_ITEM_IN_LIST_INTG1",err,error,*999)

    IF(ASSOCIATED(LIST)) THEN
      IF(list%LIST_FINISHED) THEN
        IF(list%DATA_TYPE==LIST_INTG_TYPE) THEN
!!TODO: Could search better but requires list to be sorted.
          IF(list%DATA_DIMENSION==1) THEN
            CALL LIST_SEARCH_LINEAR(list%LIST_INTG(1:list%NUMBER_IN_LIST),ITEM,LIST_ITEM,err,error,*999)
          ELSE
            CALL LIST_SEARCH_LINEAR(list%LIST_INTG2(list%KEY_DIMENSION,1:list%NUMBER_IN_LIST),ITEM,LIST_ITEM,err,error,*999)
          ENDIF
        ELSE
          localError="The list data type of "//TRIM(NumberToVString(list%DATA_TYPE,"*",err,error))// &
            & " does not match the integer type of the supplied list item."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("List has not been finished.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("List is not associated.",err,error,*999)
    ENDIF

    Call Exits("LIST_ITEM_IN_LIST_INTG1")
    RETURN
999 CALL Errors("LIST_ITEM_IN_LIST_INTG1",err,error)
    Call Exits("LIST_ITEM_IN_LIST_INTG1")
    RETURN 1
  END SUBROUTINE LIST_ITEM_IN_LIST_INTG1
  
  !
  !================================================================================================================================
  !

  !>Determines if ITEM is in the given integer LIST. If it is LIST_ITEM is the index in the list. If not LIST_ITEM is 0.
  SUBROUTINE LIST_ITEM_IN_LIST_INTG2(LIST,ITEM,LIST_ITEM,err,error,*)

    !Argument Variables
    TYPE(LIST_TYPE), POINTER, INTENT(IN) :: LIST !<The pointer to the list
    INTEGER(INTG), INTENT(IN) :: ITEM(:) !<The item to find.
    INTEGER(INTG), INTENT(OUT) :: LIST_ITEM !<On exit, the position of the item in the list. If the item does not exist then the value of 0 is returned.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    CALL Enters("LIST_ITEM_IN_LIST_INTG2",err,error,*999)

    IF(ASSOCIATED(LIST)) THEN
      IF(list%LIST_FINISHED) THEN
        IF(list%DATA_TYPE==LIST_INTG_TYPE) THEN
!!TODO: Could search better but requires list to be sorted.
          IF(list%DATA_DIMENSION==1) THEN
            CALL LIST_SEARCH_LINEAR(list%LIST_INTG(1:list%NUMBER_IN_LIST),ITEM(list%KEY_DIMENSION),LIST_ITEM,err,error,*999)
          ELSE
            CALL LIST_SEARCH_LINEAR(list%LIST_INTG2(list%KEY_DIMENSION,1:list%NUMBER_IN_LIST),ITEM(list%KEY_DIMENSION), &
              & LIST_ITEM,err,error,*999)
          ENDIF
        ELSE
          localError="The list data type of "//TRIM(NumberToVString(list%DATA_TYPE,"*",err,error))// &
            & " does not match the integer type of the supplied list item."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("List has not been finished.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("List is not associated.",err,error,*999)
    ENDIF

    Call Exits("LIST_ITEM_IN_LIST_INTG2")
    RETURN
999 CALL Errors("LIST_ITEM_IN_LIST_INTG2",err,error)
    Call Exits("LIST_ITEM_IN_LIST_INTG2")
    RETURN 1
  END SUBROUTINE LIST_ITEM_IN_LIST_INTG2
  
  !
  !================================================================================================================================
  !

  !> Determines if ITEM is in the given single precision real LIST. If it is LIST_ITEM is the index in the list. If not
  !> LIST_ITEM is 0.
  SUBROUTINE LIST_ITEM_IN_LIST_SP1(LIST,ITEM,LIST_ITEM,err,error,*)

    !Argument Variables    
    TYPE(LIST_TYPE), POINTER, INTENT(IN) :: LIST !<The pointer to the list
    REAL(SP), INTENT(IN) :: ITEM !<The item to find.
    INTEGER(INTG), INTENT(OUT) :: LIST_ITEM !<On exit, the position of the item in the list. If the item does not exist then the value of 0 is returned.     
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.    
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    CALL Enters("LIST_ITEM_IN_LIST_SP1",err,error,*999)

    IF(ASSOCIATED(LIST)) THEN
      IF(list%LIST_FINISHED) THEN
        IF(list%DATA_TYPE==LIST_SP_TYPE) THEN
!!TODO: Could search better but requires list to be sorted.
          IF(list%DATA_DIMENSION==1) THEN
            CALL LIST_SEARCH_LINEAR(list%LIST_SP(1:list%NUMBER_IN_LIST),ITEM,LIST_ITEM,err,error,*999)
          ELSE
            CALL LIST_SEARCH_LINEAR(list%LIST_SP2(lIST%KEY_DIMENSION,1:list%NUMBER_IN_LIST),ITEM,LIST_ITEM,err,error,*999)
          ENDIF
        ELSE
          localError="The list data type of "//TRIM(NumberToVString(list%DATA_TYPE,"*",err,error))// &
            & " does not match the single precision type of the supplied list item."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("List has not been finished.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("List is not associated.",err,error,*999)
    ENDIF

    Call Exits("LIST_ITEM_IN_LIST_SP1")
    RETURN
999 CALL Errors("LIST_ITEM_IN_LIST_SP1",err,error)
    Call Exits("LIST_ITEM_IN_LIST_SP1")
    RETURN 1
  END SUBROUTINE LIST_ITEM_IN_LIST_SP1
  
  !
  !================================================================================================================================
  !

  !> Determines if ITEM is in the given single precision real LIST. If it is LIST_ITEM is the index in the list. If not
  !> LIST_ITEM is 0.
  SUBROUTINE LIST_ITEM_IN_LIST_SP2(LIST,ITEM,LIST_ITEM,err,error,*)

    !Argument Variables    
    TYPE(LIST_TYPE), POINTER, INTENT(IN) :: LIST !<The pointer to the list
    REAL(SP), INTENT(IN) :: ITEM(:) !<The item to find.
    INTEGER(INTG), INTENT(OUT) :: LIST_ITEM !<On exit, the position of the item in the list. If the item does not exist then the value of 0 is returned.     
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.    
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    CALL Enters("LIST_ITEM_IN_LIST_SP2",err,error,*999)

    IF(ASSOCIATED(LIST)) THEN
      IF(list%LIST_FINISHED) THEN
        IF(list%DATA_TYPE==LIST_SP_TYPE) THEN
!!TODO: Could search better but requires list to be sorted.
          IF(list%DATA_DIMENSION==1) THEN
            CALL LIST_SEARCH_LINEAR(list%LIST_SP(1:list%NUMBER_IN_LIST),ITEM(list%KEY_DIMENSION),LIST_ITEM,err,error,*999)
          ELSE
            CALL LIST_SEARCH_LINEAR(list%LIST_SP2(lIST%KEY_DIMENSION,1:list%NUMBER_IN_LIST),ITEM(list%KEY_DIMENSION), &
              & LIST_ITEM,err,error,*999)
          ENDIF
        ELSE
          localError="The list data type of "//TRIM(NumberToVString(list%DATA_TYPE,"*",err,error))// &
            & " does not match the single precision type of the supplied list item."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("List has not been finished.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("List is not associated.",err,error,*999)
    ENDIF

    Call Exits("LIST_ITEM_IN_LIST_SP2")
    RETURN
999 CALL Errors("LIST_ITEM_IN_LIST_SP2",err,error)
    Call Exits("LIST_ITEM_IN_LIST_SP2")
    RETURN 1
  END SUBROUTINE LIST_ITEM_IN_LIST_SP2
  
  !
  !================================================================================================================================
  !

  !> Determines if ITEM is in the given double precision real LIST. If it is LIST_ITEM is the index in the list. If not
  !> LIST_ITEM is 0.
  SUBROUTINE LIST_ITEM_IN_LIST_DP1(LIST,ITEM,LIST_ITEM,err,error,*)

    !Argument Variables
    TYPE(LIST_TYPE), POINTER, INTENT(IN) :: LIST !<The pointer to the list
    REAL(DP), INTENT(IN) :: ITEM  !<The item to find.
    INTEGER(INTG), INTENT(OUT) :: LIST_ITEM !<On exit, the position of the item in the list. If the item does not exist then the value of 0 is returned.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    CALL Enters("LIST_ITEM_IN_LIST_DP1",err,error,*999)

    IF(ASSOCIATED(LIST)) THEN
      IF(list%LIST_FINISHED) THEN
        IF(list%DATA_TYPE==LIST_DP_TYPE) THEN
!!TODO: Could search better but requires list to be sorted.
          IF(list%DATA_DIMENSION==1) THEN
            CALL LIST_SEARCH_LINEAR(list%LIST_DP(1:list%NUMBER_IN_LIST),ITEM,LIST_ITEM,err,error,*999)
          ELSE
            CALL LIST_SEARCH_LINEAR(list%LIST_DP2(list%KEY_DIMENSION,1:list%NUMBER_IN_LIST),ITEM,LIST_ITEM,err,error,*999)
         ENDIF
        ELSE
          localError="The list data type of "//TRIM(NumberToVString(list%DATA_TYPE,"*",err,error))// &
            & " does not match the single precision type of the supplied list item."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("List has not been finished.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("List is not associated.",err,error,*999)
    ENDIF

    Call Exits("LIST_ITEM_IN_LIST_DP1")
    RETURN
999 CALL Errors("LIST_ITEM_IN_LIST_DP1",err,error)
    Call Exits("LIST_ITEM_IN_LIST_DP1")
    RETURN 1
  END SUBROUTINE LIST_ITEM_IN_LIST_DP1

  !
  !================================================================================================================================
  !

  !> Determines if ITEM is in the given double precision real LIST. If it is LIST_ITEM is the index in the list. If not
  !> LIST_ITEM is 0.
  SUBROUTINE LIST_ITEM_IN_LIST_DP2(LIST,ITEM,LIST_ITEM,err,error,*)

    !Argument Variables
    TYPE(LIST_TYPE), POINTER, INTENT(IN) :: LIST !<The pointer to the list
    REAL(DP), INTENT(IN) :: ITEM(:)  !<The item to find.
    INTEGER(INTG), INTENT(OUT) :: LIST_ITEM !<On exit, the position of the item in the list. If the item does not exist then the value of 0 is returned.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    CALL Enters("LIST_ITEM_IN_LIST_DP2",err,error,*999)

    IF(ASSOCIATED(LIST)) THEN
      IF(list%LIST_FINISHED) THEN
        IF(list%DATA_TYPE==LIST_DP_TYPE) THEN
!!TODO: Could search better but requires list to be sorted.
          IF(list%DATA_DIMENSION==1) THEN
            CALL LIST_SEARCH_LINEAR(list%LIST_DP(1:list%NUMBER_IN_LIST),ITEM(list%KEY_DIMENSION),LIST_ITEM,err,error,*999)
          ELSE
            CALL LIST_SEARCH_LINEAR(list%LIST_DP2(list%KEY_DIMENSION,1:list%NUMBER_IN_LIST),ITEM(list%KEY_DIMENSION), &
              & LIST_ITEM,err,error,*999)
         ENDIF
        ELSE
          localError="The list data type of "//TRIM(NumberToVString(list%DATA_TYPE,"*",err,error))// &
            & " does not match the single precision type of the supplied list item."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("List has not been finished.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("List is not associated.",err,error,*999)
    ENDIF

    Call Exits("LIST_ITEM_IN_LIST_DP2")
    RETURN
999 CALL Errors("LIST_ITEM_IN_LIST_DP2",err,error)
    Call Exits("LIST_ITEM_IN_LIST_DP2")
    RETURN 1
  END SUBROUTINE LIST_ITEM_IN_LIST_DP2

  !
  !================================================================================================================================
  !
  
  !>Deletes the item given by the listItem index from the given list.
  SUBROUTINE List_ItemDelete(list,listItem,err,error,*)

    !Argument Variables
    TYPE(LIST_TYPE), POINTER, INTENT(IN) :: list !<The pointer to the list
    INTEGER(INTG), INTENT(IN) :: listItem !<The position in the list to delete.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    CALL Enters("List_ItemDelete",err,error,*999)

    IF(ASSOCIATED(list)) THEN
      IF(list%LIST_FINISHED) THEN
        IF(listItem>=1.AND.listItem<=list%NUMBER_IN_LIST) THEN
          IF(list%DATA_DIMENSION==1) THEN
            SELECT CASE(list%DATA_TYPE)
            CASE(LIST_INTG_TYPE)
              list%LIST_INTG(1:listItem-1)=list%LIST_INTG(1:listItem-1)
              list%LIST_INTG(listItem:list%NUMBER_IN_LIST-1)=list%LIST_INTG(listItem+1:list%NUMBER_IN_LIST)
            CASE(LIST_SP_TYPE)
              list%LIST_SP(1:listItem-1)=list%LIST_SP(1:listItem-1)
              list%LIST_SP(listItem:list%NUMBER_IN_LIST-1)=list%LIST_SP(listItem+1:list%NUMBER_IN_LIST)
            CASE(LIST_DP_TYPE)
              list%LIST_DP(1:listItem-1)=list%LIST_DP(1:listItem-1)
              list%LIST_DP(listItem:list%NUMBER_IN_LIST-1)=list%LIST_DP(listItem+1:list%NUMBER_IN_LIST)
            CASE DEFAULT
              localError="The list data type of "//TRIM(NumberToVString(list%DATA_TYPE,"*",err,error))//" is invalid."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          ELSE
            SELECT CASE(list%DATA_TYPE)
            CASE(LIST_INTG_TYPE)
              list%LIST_INTG2(:,1:listItem-1)=list%LIST_INTG2(:,1:listItem-1)
              list%LIST_INTG2(:,listItem:list%NUMBER_IN_LIST-1)=list%LIST_INTG2(:,listItem+1:list%NUMBER_IN_LIST)
            CASE(LIST_SP_TYPE)
              list%LIST_SP2(:,1:listItem-1)=list%LIST_SP2(:,1:listItem-1)
              list%LIST_SP2(:,listItem:list%NUMBER_IN_LIST-1)=list%LIST_SP2(:,listItem+1:list%NUMBER_IN_LIST)
            CASE(LIST_DP_TYPE)
              list%LIST_DP2(:,1:listItem-1)=list%LIST_DP2(:,1:listItem-1)
              list%LIST_DP2(:,listItem:list%NUMBER_IN_LIST-1)=list%LIST_DP2(:,listItem+1:list%NUMBER_IN_LIST)
            CASE DEFAULT
              localError="The list data type of "//TRIM(NumberToVString(list%DATA_TYPE,"*",err,error))//" is invalid."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          ENDIF
          list%NUMBER_IN_LIST=list%NUMBER_IN_LIST-1
        ELSE
          localError="The specified list item of "//TRIM(NumberToVString(listItem,"*",err,error))// &
            & " is invalid. The item must be >= 1 and <= "//TRIM(NumberToVString(list%NUMBER_IN_LIST,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("List has not been finished.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("List is not associated.",err,error,*999)
    ENDIF

    Call Exits("List_ItemDelete")
    RETURN
999 CALL Errors("List_ItemDelete",err,error)
    Call Exits("List_ItemDelete")
    RETURN 1
  END SUBROUTINE List_ItemDelete
  
 !
  !================================================================================================================================
  !
  
  !>Deletes the item given by the LIST_ITEM index from the given list.
  SUBROUTINE LIST_ITEM_DELETE(LIST,LIST_ITEM,err,error,*)

    !Argument Variables
    TYPE(LIST_TYPE), POINTER, INTENT(IN) :: LIST !<The pointer to the list
    INTEGER(INTG), INTENT(IN) :: LIST_ITEM !<The position in the list to delete.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL List_ItemDelete(LIST,LIST_ITEM,err,error,*999)
    
    RETURN
999 RETURN 1
  END SUBROUTINE LIST_ITEM_DELETE
  
  !
  !================================================================================================================================
  !

  !>Sets/changes the key dimension (i.e., the dimension for searching and sorting) for a list
  SUBROUTINE List_KeyDimensionSet(list,keyDimension,err,error,*)

    !Argument Variables
    TYPE(LIST_TYPE), POINTER, INTENT(IN) :: list !<A pointer to the list
    INTEGER(INTG), INTENT(IN) :: keyDimension !<The key dimension to set. Must be greater than zero and <= the data dimension.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    CALL Enters("List_KeyDimensionSet",err,error,*999)

    IF(ASSOCIATED(list)) THEN
      IF(keyDimension>0.AND.keyDimension<=list%DATA_DIMENSION) THEN
        list%KEY_DIMENSION=keyDimension
      ELSE
        localError="The specified key dimension of "//TRIM(NumberToVString(keyDimension,"*",err,error))// &
          & " is invalid. The key dimension must be > 0 and <= "// &
          & TRIM(NumberToVString(list%DATA_DIMENSION,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("List is not associated.",err,error,*999)
    ENDIF

    Call Exits("List_KeyDimensionSet")
    RETURN
999 CALL Errors("List_KeyDimensionSet",err,error)
    Call Exits("List_KeyDimensionSet")
    RETURN 1
  END SUBROUTINE List_KeyDimensionSet

  !
  !================================================================================================================================
  !

  !>Sets/changes the key dimension (i.e., the dimension for searching and sorting) for a list
  SUBROUTINE LIST_KEY_DIMENSION_SET(LIST,KEY_DIMENSION,err,error,*)

    !Argument Variables
    TYPE(LIST_TYPE), POINTER, INTENT(IN) :: LIST !<A pointer to the list
    INTEGER(INTG), INTENT(IN) :: KEY_DIMENSION !<The key dimension to set. Must be greater than zero and <= the data dimension.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL List_KeyDimensionSet(LIST,KEY_DIMENSION,err,error,*999)
    
    RETURN
999 RETURN 1
  END SUBROUTINE LIST_KEY_DIMENSION_SET

  !
  !================================================================================================================================
  !

  !>Gets the current number of items in a list
  SUBROUTINE List_NumberOfItemsGet(list,numberOfItems,err,error,*)
      
    !Argument variables
    TYPE(LIST_TYPE), POINTER, INTENT(IN) :: list !<A pointer to the list 
    INTEGER(INTG), INTENT(OUT) :: numberOfItems !<On exit, the current number of items in the list
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    
    CALL Enters("List_NumberOfItemsGet",err,error,*999)

    IF(ASSOCIATED(list)) THEN
      IF(list%LIST_FINISHED) THEN
        numberOfItems=list%NUMBER_IN_LIST
      ELSE
        CALL FlagError("List has not been finished.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("List is not associated.",err,error,*999)
    ENDIF
    
    Call Exits("List_NumberOfItemsGet")
    RETURN
999 CALL Errors("List_NumberOfItemsGet",err,error)
    Call Exits("List_NumberOfItemsGet")
    RETURN 1
  END SUBROUTINE List_NumberOfItemsGet
  
  !
  !================================================================================================================================
  !

  !>Gets the current number of items in a list
  SUBROUTINE LIST_NUMBER_OF_ITEMS_GET(LIST,NUMBER_OF_ITEMS,err,error,*)
      
    !Argument variables
    TYPE(LIST_TYPE), POINTER, INTENT(IN) :: LIST !<A pointer to the list 
    INTEGER(INTG), INTENT(OUT) :: NUMBER_OF_ITEMS !<On exit, the current number of items in the list
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables

    CALL List_NumberOfItemsGet(LIST,NUMBER_OF_ITEMS,err,error,*999)
    
    RETURN
999 RETURN 1
  END SUBROUTINE LIST_NUMBER_OF_ITEMS_GET
  
  !
  !================================================================================================================================
  !

  !>Detaches the list values from an integer list of data dimension 1 and returns them as an array of base type
  !>before destroying the list. The LIST_VALUES array must not be allocated on entry. It is up to the user to then deallocate
  !>the returned list memory.
  SUBROUTINE LIST_DETACH_AND_DESTROY_INTG1(LIST,NUMBER_IN_LIST,LIST_VALUES,err,error,*)

    !Argument Variables
    TYPE(LIST_TYPE), POINTER, INTENT(INOUT) :: LIST !<The pointer to the list
    INTEGER(INTG), INTENT(OUT) :: NUMBER_IN_LIST !<On exit, the number in the list that has been detached.
    INTEGER(INTG), ALLOCATABLE, INTENT(INOUT) :: LIST_VALUES(:) !<On exit, the detached list. Must not be allocated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    CALL Enters("LIST_DETACH_AND_DESTROY_INTG1",err,error,*999)

    IF(ASSOCIATED(LIST)) THEN
      IF(list%LIST_FINISHED) THEN
        IF(ALLOCATED(LIST_VALUES)) THEN
          CALL FlagError("List values is allocated.",err,error,*999)
        ELSE
          IF(list%DATA_TYPE==LIST_INTG_TYPE) THEN
            IF(list%DATA_DIMENSION==1) THEN
              NUMBER_IN_LIST=list%NUMBER_IN_LIST
              !Note this will return more memory as the list will be bigger. Maybe copy to an array the correct size?
              CALL MOVE_ALLOC(list%LIST_INTG,LIST_VALUES)
              CALL LIST_FINALISE(LIST,err,error,*999)
            ELSE
              localError="Invalid data dimension. The supplied data dimension is 1 and the list data dimension is "// &
                & TRIM(NumberToVString(list%DATA_DIMENSION,"*",err,error))//"."
              CALL FlagError(localError,err,error,*999)
            ENDIF
          ELSE
            localError="The list data type of "//TRIM(NumberToVString(list%DATA_TYPE,"*",err,error))// &
              & " does not match the integer type of the supplied list values item."
            CALL FlagError(localError,err,error,*999)
          ENDIF
        ENDIF
      ELSE
        CALL FlagError("List has not been finished.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("List is not associated.",err,error,*999)
    ENDIF
    
    Call Exits("LIST_DETACH_AND_DESTROY_INTG1")
    RETURN
999 CALL Errors("LIST_DETACH_AND_DESTROY_INTG1",err,error)
    Call Exits("LIST_DETACH_AND_DESTROY_INTG1")
    RETURN 1
  END SUBROUTINE LIST_DETACH_AND_DESTROY_INTG1

  !
  !================================================================================================================================
  !

  !>Detaches the list values from an integer list of data dimension > 1 and returns them as an array of base type
  !>before destroying the list. The LIST_VALUES array must not be allocated on entry. It is up to the user to then deallocate
  !>the returned list memory.
  SUBROUTINE LIST_DETACH_AND_DESTROY_INTG2(LIST,NUMBER_IN_LIST,LIST_VALUES,err,error,*)

    !Argument Variables
    TYPE(LIST_TYPE), POINTER, INTENT(INOUT) :: LIST !<The pointer to the list
    INTEGER(INTG), INTENT(OUT) :: NUMBER_IN_LIST !<On exit, the number in the list that has been detached.
    INTEGER(INTG), ALLOCATABLE, INTENT(INOUT) :: LIST_VALUES(:,:) !<On exit, the detached list. Must not be allocated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    CALL Enters("LIST_DETACH_AND_DESTROY_INTG2",err,error,*999)

    IF(ASSOCIATED(LIST)) THEN
      IF(list%LIST_FINISHED) THEN
        IF(ALLOCATED(LIST_VALUES)) THEN
          CALL FlagError("List values is allocated.",err,error,*999)
        ELSE
          IF(list%DATA_TYPE==LIST_INTG_TYPE) THEN
            IF(list%DATA_DIMENSION>1) THEN
              NUMBER_IN_LIST=list%NUMBER_IN_LIST
              !Note this will return more memory as the list will be bigger. Maybe copy to an array the correct size?
              CALL MOVE_ALLOC(list%LIST_INTG2,LIST_VALUES)
              CALL LIST_FINALISE(LIST,err,error,*999)
            ELSE
              CALL FlagError("Invalid data dimension. The supplied data dimension is > 1 and the list data dimension is 1.", &
                & err,error,*999)
            ENDIF
          ELSE
            localError="The list data type of "//TRIM(NumberToVString(list%DATA_TYPE,"*",err,error))// &
              & " does not match the integer type of the supplied list values item."
            CALL FlagError(localError,err,error,*999)
          ENDIF
        ENDIF
      ELSE
        CALL FlagError("List has not been finished.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("List is not associated.",err,error,*999)
    ENDIF
    
    Call Exits("LIST_DETACH_AND_DESTROY_INTG2")
    RETURN
999 CALL Errors("LIST_DETACH_AND_DESTROY_INTG2",err,error)
    Call Exits("LIST_DETACH_AND_DESTROY_INTG2")
    RETURN 1
  END SUBROUTINE LIST_DETACH_AND_DESTROY_INTG2

  !
  !================================================================================================================================
  !

  !>Detaches the list values from a single precision real list of data dimension 1 and returns them as an array
  !>of base type before destroying the list. The LIST_VALUES array must not be allocated on entry. It is up to the user to
  !>then deallocate the returned list memory.
  SUBROUTINE LIST_DETACH_AND_DESTROY_SP1(LIST,NUMBER_IN_LIST,LIST_VALUES,err,error,*)

    !Argument Variables
    TYPE(LIST_TYPE), POINTER, INTENT(INOUT) :: LIST !<The pointer to the list
    INTEGER(INTG), INTENT(OUT) :: NUMBER_IN_LIST !<On exit, the number in the list that has been detached.
    REAL(SP), ALLOCATABLE, INTENT(INOUT) :: LIST_VALUES(:) !<On exit, the detached list. Must not be allocated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    CALL Enters("LIST_DETACH_AND_DESTROY_SP1",err,error,*999)

    IF(ASSOCIATED(LIST)) THEN
      IF(list%LIST_FINISHED) THEN
        IF(ALLOCATED(LIST_VALUES)) THEN
          CALL FlagError("List values is associated.",err,error,*999)
        ELSE
          IF(list%DATA_TYPE==LIST_SP_TYPE) THEN
            IF(list%DATA_DIMENSION==1) THEN
              NUMBER_IN_LIST=list%NUMBER_IN_LIST
              !Note this will return more memory as the list will be bigger. Maybe copy to an array the correct size?
              CALL MOVE_ALLOC(list%LIST_SP,LIST_VALUES)
              CALL LIST_FINALISE(LIST,err,error,*999)
            ELSE
              localError="Invalid data dimension. The supplied data dimension is 1 and the list data dimension is "// &
                & TRIM(NumberToVString(list%DATA_DIMENSION,"*",err,error))//"."
              CALL FlagError(localError,err,error,*999)
            ENDIF
          ELSE
            localError="The list data type of "//TRIM(NumberToVString(list%DATA_TYPE,"*",err,error))// &
              & " does not match the single precision type of the supplied list values item."
            CALL FlagError(localError,err,error,*999)
          ENDIF
        ENDIF
      ELSE
        CALL FlagError("List has not been finished.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("List is not associated.",err,error,*999)
    ENDIF
    
    Call Exits("LIST_DETACH_AND_DESTROY_SP1")
    RETURN
999 CALL Errors("LIST_DETACH_AND_DESTROY_SP1",err,error)
    Call Exits("LIST_DETACH_AND_DESTROY_SP1")
    RETURN 1
  END SUBROUTINE LIST_DETACH_AND_DESTROY_SP1
  !
  !================================================================================================================================
  !
  
  !>Detaches the list values from a single precision real list of data dimension > 1 and returns them as an array
  !>of base type before destroying the list. The LIST_VALUES array must not be allocated on entry. It is up to the user to
  !>then deallocate the returned list memory.
  SUBROUTINE LIST_DETACH_AND_DESTROY_SP2(LIST,NUMBER_IN_LIST,LIST_VALUES,err,error,*)

    !Argument Variables
    TYPE(LIST_TYPE), POINTER, INTENT(INOUT) :: LIST !<The pointer to the list
    INTEGER(INTG), INTENT(OUT) :: NUMBER_IN_LIST !<On exit, the number in the list that has been detached.
    REAL(SP), ALLOCATABLE, INTENT(INOUT) :: LIST_VALUES(:,:) !<On exit, the detached list. Must not be allocated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    CALL Enters("LIST_DETACH_AND_DESTROY_SP2",err,error,*999)

    IF(ASSOCIATED(LIST)) THEN
      IF(list%LIST_FINISHED) THEN
        IF(ALLOCATED(LIST_VALUES)) THEN
          CALL FlagError("List values is associated.",err,error,*999)
        ELSE
          IF(list%DATA_TYPE==LIST_SP_TYPE) THEN
            IF(list%DATA_DIMENSION>1) THEN
              NUMBER_IN_LIST=list%NUMBER_IN_LIST
              !Note this will return more memory as the list will be bigger. Maybe copy to an array the correct size?
              CALL MOVE_ALLOC(list%LIST_SP2,LIST_VALUES)
              CALL LIST_FINALISE(LIST,err,error,*999)
            ELSE
              CALL FlagError("Invalid data dimension. The supplied data dimension is > 1 and the list data dimension is 1.", &
                & err,error,*999)
            ENDIF
          ELSE
            localError="The list data type of "//TRIM(NumberToVString(list%DATA_TYPE,"*",err,error))// &
              & " does not match the single precision type of the supplied list values item."
            CALL FlagError(localError,err,error,*999)
          ENDIF
        ENDIF
      ELSE
        CALL FlagError("List has not been finished.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("List is not associated.",err,error,*999)
    ENDIF
    
    Call Exits("LIST_DETACH_AND_DESTROY_SP2")
    RETURN
999 CALL Errors("LIST_DETACH_AND_DESTROY_SP2",err,error)
    Call Exits("LIST_DETACH_AND_DESTROY_SP2")
    RETURN 1
  END SUBROUTINE LIST_DETACH_AND_DESTROY_SP2

  !
  !================================================================================================================================
  !
  
  !>Detaches the list values from a double precision real list of data dimension 1 and returns them as an array
  !>of base type before destroying the list. The LIST_VALUES array must not be allocated on entry. It is up to the user
  !>to then deallocate the returned list memory.
  SUBROUTINE LIST_DETACH_AND_DESTROY_DP1(LIST,NUMBER_IN_LIST,LIST_VALUES,err,error,*)

    !Argument Variables
    TYPE(LIST_TYPE), POINTER, INTENT(INOUT) :: LIST !<The pointer to the list
    INTEGER(INTG), INTENT(OUT) :: NUMBER_IN_LIST !<On exit, the number in the list that has been detached.
    REAL(DP), ALLOCATABLE, INTENT(INOUT) :: LIST_VALUES(:) !<On exit, the detached list. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    CALL Enters("LIST_DETACH_AND_DESTROY_DP1",err,error,*999)

    IF(ASSOCIATED(LIST)) THEN
      IF(list%LIST_FINISHED) THEN
        IF(ALLOCATED(LIST_VALUES)) THEN
          CALL FlagError("List values is associated.",err,error,*999)
        ELSE
          IF(list%DATA_TYPE==LIST_DP_TYPE) THEN
            IF(list%DATA_DIMENSION==1) THEN
              NUMBER_IN_LIST=list%NUMBER_IN_LIST
              !Note this will return more memory as the list will be bigger. Maybe copy to an array the correct size?
              CALL MOVE_ALLOC(list%LIST_DP,LIST_VALUES)
              CALL LIST_FINALISE(LIST,err,error,*999)
            ELSE
              localError="Invalid data dimension. The supplied data dimension is 1 and the list data dimension is "// &
                & TRIM(NumberToVString(list%DATA_DIMENSION,"*",err,error))//"."
              CALL FlagError(localError,err,error,*999)
            ENDIF
          ELSE
            localError="The list data type of "//TRIM(NumberToVString(list%DATA_TYPE,"*",err,error))// &
              & " does not match the double precision type of the supplied list values item."
            CALL FlagError(localError,err,error,*999)
          ENDIF
        ENDIF
      ELSE
        CALL FlagError("List has not been finished.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("List is not associated.",err,error,*999)
    ENDIF
    
    Call Exits("LIST_DETACH_AND_DESTROY_DP1")
    RETURN
999 CALL Errors("LIST_DETACH_AND_DESTROY_DP1",err,error)
    Call Exits("LIST_DETACH_AND_DESTROY_DP1")
    RETURN 1
  END SUBROUTINE LIST_DETACH_AND_DESTROY_DP1

  !
  !================================================================================================================================
  !
  
  !>Detaches the list values from a double precision real list of data dimension > 1 and returns them as an array
  !>of base type before destroying the list. The LIST_VALUES array must not be allocated on entry. It is up to the user
  !>to then deallocate the returned list memory.
  SUBROUTINE LIST_DETACH_AND_DESTROY_DP2(LIST,NUMBER_IN_LIST,LIST_VALUES,err,error,*)

    !Argument Variables
    TYPE(LIST_TYPE), POINTER, INTENT(INOUT) :: LIST !<The pointer to the list
    INTEGER(INTG), INTENT(OUT) :: NUMBER_IN_LIST !<On exit, the number in the list that has been detached.
    REAL(DP), ALLOCATABLE, INTENT(INOUT) :: LIST_VALUES(:,:) !<On exit, the detached list. Must not be allocated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    CALL Enters("LIST_DETACH_AND_DESTROY_DP2",err,error,*999)

    IF(ASSOCIATED(LIST)) THEN
      IF(list%LIST_FINISHED) THEN
        IF(ALLOCATED(LIST_VALUES)) THEN
          CALL FlagError("List values is associated.",err,error,*999)
        ELSE
          IF(list%DATA_TYPE==LIST_DP_TYPE) THEN
            IF(list%DATA_DIMENSION>1) THEN
              NUMBER_IN_LIST=list%NUMBER_IN_LIST
              !Note this will return more memory as the list will be bigger. Maybe copy to an array the correct size?
              CALL MOVE_ALLOC(list%LIST_DP2,LIST_VALUES)
              CALL LIST_FINALISE(LIST,err,error,*999)
            ELSE
              CALL FlagError("Invalid data dimension. The supplied data dimension is > 1 and the list data dimension is 1.", &
                & err,error,*999)
            ENDIF
          ELSE
            localError="The list data type of "//TRIM(NumberToVString(list%DATA_TYPE,"*",err,error))// &
              & " does not match the double precision type of the supplied list values item."
            CALL FlagError(localError,err,error,*999)
          ENDIF
        ENDIF
      ELSE
        CALL FlagError("List has not been finished.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("List is not associated.",err,error,*999)
    ENDIF
    
    Call Exits("LIST_DETACH_AND_DESTROY_DP2")
    RETURN
999 CALL Errors("LIST_DETACH_AND_DESTROY_DP2",err,error)
    Call Exits("LIST_DETACH_AND_DESTROY_DP2")
    RETURN 1
  END SUBROUTINE LIST_DETACH_AND_DESTROY_DP2

  !
  !================================================================================================================================
  !

  !>Removes duplicate entries from a list. A side effect of this is that the list is sorted.
  SUBROUTINE List_RemoveDuplicates(list,err,error,*)

    !Argument Variables
    TYPE(LIST_TYPE), POINTER, INTENT(INOUT) :: list !<The pointer to the list
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: i,j,numberRemoved
    LOGICAL :: sameValue
    TYPE(VARYING_STRING) :: localError

    CALL Enters("List_RemoveDuplicates",err,error,*999)

    IF(ASSOCIATED(list)) THEN
      IF(list%LIST_FINISHED) THEN
        IF(list%NUMBER_IN_LIST>0) THEN
          IF(list%DATA_DIMENSION==1) THEN
            SELECT CASE(list%DATA_TYPE)
            CASE(LIST_INTG_TYPE)              
              CALL List_Sort(list%LIST_INTG(1:list%NUMBER_IN_LIST),err,error,*999)
              i=1
              DO WHILE(i<=list%NUMBER_IN_LIST)
                !Find the extent of duplicate values if any
                j=i+1
                sameValue=.TRUE.
                DO WHILE(j<=list%NUMBER_IN_LIST.AND.sameValue)
                  IF(list%LIST_INTG(j)==list%LIST_INTG(i)) THEN
                    j=j+1
                  ELSE
                    sameValue=.FALSE.
                  ENDIF
                ENDDO !j
                IF(j>i+1.OR.sameValue) THEN
                  !We have duplicates so remove them
                  IF(sameValue) THEN
                    !Duplicates to the end of the list so just set the number in the list
                    list%NUMBER_IN_LIST=i
                  ELSE
                    numberRemoved=j-i-1
                    list%LIST_INTG(i+1:list%NUMBER_IN_LIST-numberRemoved)=list%LIST_INTG(j:list%NUMBER_IN_LIST)
                    list%NUMBER_IN_LIST=list%NUMBER_IN_LIST-numberRemoved
                  ENDIF
                ENDIF
                i=i+1
              ENDDO !i
            CASE(LIST_SP_TYPE)
              CALL LIST_SORT(list%LIST_SP(1:list%NUMBER_IN_LIST),err,error,*999)
              i=1
              DO WHILE(i<=list%NUMBER_IN_LIST)
                !Find the extent of duplicate values if any
                j=i+1
                sameValue=.TRUE.
                DO WHILE(j<=list%NUMBER_IN_LIST.AND.sameValue)
                  IF(ABS(list%LIST_SP(j)-list%LIST_SP(i))<=ZERO_TOLERANCE_SP) THEN
                    j=j+1
                  ELSE
                    sameValue=.FALSE.
                  ENDIF
                ENDDO !j
                IF(j>i+1.OR.sameValue) THEN
                  !We have duplicates so remove them
                  IF(sameValue) THEN
                    !Duplicates to the end of the list so just set the number in the list
                    list%NUMBER_IN_LIST=i
                  ELSE
                    numberRemoved=j-i-1
                    list%LIST_SP(i+1:list%NUMBER_IN_LIST-numberRemoved)=list%LIST_SP(j:list%NUMBER_IN_LIST)
                    list%NUMBER_IN_LIST=list%NUMBER_IN_LIST-numberRemoved
                  ENDIF
                ENDIF
                i=i+1
              ENDDO !i
            CASE(LIST_DP_TYPE)
              CALL LIST_SORT(list%LIST_DP(1:list%NUMBER_IN_LIST),err,error,*999)
              i=1
              DO WHILE(i<=list%NUMBER_IN_LIST)
                !Find the extent of duplicate values if any
                j=i+1
                sameValue=.TRUE.
                DO WHILE(j<=list%NUMBER_IN_LIST.AND.sameValue)
                  IF(ABS(list%LIST_DP(j)-list%LIST_DP(i))<=ZERO_TOLERANCE) THEN
                    j=j+1
                  ELSE
                    sameValue=.FALSE.
                  ENDIF
                ENDDO !j
                IF(j>i+1.OR.sameValue) THEN
                  !We have duplicates so remove them
                  IF(sameValue) THEN
                    !Duplicates to the end of the list so just set the number in the list
                    list%NUMBER_IN_LIST=i
                  ELSE
                    numberRemoved=j-i-1
                    list%LIST_DP(i+1:list%NUMBER_IN_LIST-numberRemoved)=list%LIST_DP(j:list%NUMBER_IN_LIST)
                    list%NUMBER_IN_LIST=list%NUMBER_IN_LIST-numberRemoved
                  ENDIF
                ENDIF
                i=i+1
              ENDDO !i
            CASE DEFAULT
              localError="The list data type of "//TRIM(NumberToVString(list%DATA_TYPE,"*",err,error))//" is invalid."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          ELSE
            SELECT CASE(list%DATA_TYPE)
            CASE(LIST_INTG_TYPE)              
              CALL LIST_SORT(list%LIST_INTG2(:,1:list%NUMBER_IN_LIST),list%KEY_DIMENSION,err,error,*999)
              i=1
              DO WHILE(i<=list%NUMBER_IN_LIST)
                !Find the extent of duplicate values if any
                j=i+1
                sameValue=.TRUE.
                DO WHILE(j<=list%NUMBER_IN_LIST.AND.sameValue)
                  IF(list%LIST_INTG2(list%KEY_DIMENSION,j)==list%LIST_INTG2(list%KEY_DIMENSION,i)) THEN
                    j=j+1
                  ELSE
                    sameValue=.FALSE.
                  ENDIF
                ENDDO !j
                IF(j>i+1.OR.sameValue) THEN
                  !We have duplicates so remove them
                  IF(sameValue) THEN
                    !Duplicates to the end of the list so just set the number in the list
                    list%NUMBER_IN_LIST=i
                  ELSE
                    numberRemoved=j-i-1
                    list%LIST_INTG2(:,i+1:list%NUMBER_IN_LIST-numberRemoved)=list%LIST_INTG2(:,j:list%NUMBER_IN_LIST)
                    list%NUMBER_IN_LIST=list%NUMBER_IN_LIST-numberRemoved
                  ENDIF
                ENDIF
                i=i+1
              ENDDO !i
            CASE(LIST_SP_TYPE)
              CALL LIST_SORT(list%LIST_SP2(:,1:list%NUMBER_IN_LIST),list%KEY_DIMENSION,err,error,*999)
              i=1
              DO WHILE(i<=list%NUMBER_IN_LIST)
                !Find the extent of duplicate values if any
                j=i+1
                sameValue=.TRUE.
                DO WHILE(j<=list%NUMBER_IN_LIST.AND.sameValue)
                  IF(ABS(list%LIST_SP2(list%KEY_DIMENSION,j)-list%LIST_SP2(list%KEY_DIMENSION,i))<=ZERO_TOLERANCE_SP) THEN
                    j=j+1
                  ELSE
                    sameValue=.FALSE.
                  ENDIF
                ENDDO !j
                IF(j>i+1.OR.sameValue) THEN
                  !We have duplicates so remove them
                  IF(sameValue) THEN
                    !Duplicates to the end of the list so just set the number in the list
                    list%NUMBER_IN_LIST=i
                  ELSE
                    numberRemoved=j-i-1
                    list%LIST_SP2(:,i+1:list%NUMBER_IN_LIST-numberRemoved)=list%LIST_SP2(:,j:list%NUMBER_IN_LIST)
                    list%NUMBER_IN_LIST=list%NUMBER_IN_LIST-numberRemoved
                  ENDIF
                ENDIF
                i=i+1
              ENDDO !i
            CASE(LIST_DP_TYPE)
              CALL LIST_SORT(list%LIST_DP2(:,1:list%NUMBER_IN_LIST),list%KEY_DIMENSION,err,error,*999)
              i=1
              DO WHILE(i<=list%NUMBER_IN_LIST)
                !Find the extent of duplicate values if any
                j=i+1
                sameValue=.TRUE.
                DO WHILE(j<=list%NUMBER_IN_LIST.AND.sameValue)
                  IF(ABS(list%LIST_DP2(list%KEY_DIMENSION,j)-list%LIST_DP2(list%KEY_DIMENSION,i))<=ZERO_TOLERANCE) THEN
                    j=j+1
                  ELSE
                    sameValue=.FALSE.
                  ENDIF
                ENDDO !j
                IF(j>i+1.OR.sameValue) THEN
                  !We have duplicates so remove them
                  IF(sameValue) THEN
                    !Duplicates to the end of the list so just set the number in the list
                    list%NUMBER_IN_LIST=i
                  ELSE
                    numberRemoved=j-i-1
                    list%LIST_DP2(:,i+1:list%NUMBER_IN_LIST-numberRemoved)=list%LIST_DP2(:,j:list%NUMBER_IN_LIST)
                    list%NUMBER_IN_LIST=list%NUMBER_IN_LIST-numberRemoved
                  ENDIF
                ENDIF
                i=i+1
              ENDDO !i
            CASE DEFAULT
              localError="The list data type of "//TRIM(NumberToVString(list%DATA_TYPE,"*",err,error))//" is invalid."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          ENDIF
        ENDIF
      ELSE
        CALL FlagError("List has not been finished.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("List is not associated.",err,error,*999)
    ENDIF
  
    Call Exits("List_RemoveDuplicates")
    RETURN
999 CALL Errors("List_RemoveDuplicates",err,error)
    Call Exits("List_RemoveDuplicates")
    RETURN 1
  END SUBROUTINE List_RemoveDuplicates

 !
  !================================================================================================================================
  !

  !>Removes duplicate entries from a list. A side effect of this is that the list is sorted.
  SUBROUTINE LIST_REMOVE_DUPLICATES(LIST,err,error,*)

    !Argument Variables
    TYPE(LIST_TYPE), POINTER, INTENT(INOUT) :: LIST !<The pointer to the list
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL List_RemoveDuplicates(list,err,error,*999)
    
    RETURN
999 RETURN 1
  END SUBROUTINE LIST_REMOVE_DUPLICATES

  !
  !================================================================================================================================
  !

  !>Searches an integer array list A for VALUE. If the search is successful POSITION contains the index of the position of VALUE in the list otherwise POSITION is zero.
  SUBROUTINE LIST_SEARCH_INTG_ARRAY(A,VALUE,POSITION,err,error,*)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: A(:) !<The list to search
    INTEGER(INTG), INTENT(IN) :: VALUE !<The value to search for
    INTEGER(INTG), INTENT(OUT) :: POSITION !<On exit, the position of value in the list. If value does not exist in the list the returned position is zero. 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    
    CALL Enters("LIST_SEARCH_INTG_ARRAY",err,error,*999)

    !Default search method is a linear search
    CALL LIST_SEARCH_LINEAR(A,VALUE,POSITION,err,error,*999)
    
    Call Exits("LIST_SEARCH_INTG_ARRAY")
    RETURN
999 CALL Errors("LIST_SEARCH_INTG_ARRAY",err,error)
    Call Exits("LIST_SEARCH_INTG_ARRAY")
    RETURN 1
  END SUBROUTINE LIST_SEARCH_INTG_ARRAY
  
  !
  !================================================================================================================================
  !

  !>Searches an integer array list A for VALUE. If the search is successful POSITION contains the index of the position of VALUE in the list otherwise POSITION is zero.
  SUBROUTINE LIST_SEARCH_C_INT_ARRAY(A,VALUE,POSITION,err,error,*)
  
    !Argument variables
    INTEGER(C_INT), INTENT(IN) :: A(:) !<The list to search
    INTEGER(C_INT), INTENT(IN) :: VALUE !<The value to search for
    INTEGER(INTG), INTENT(OUT) :: POSITION !<On exit, the position of value in the list. If value does not exist in the list the returned position is zero. 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    
    CALL Enters("LIST_SEARCH_C_INT_ARRAY",err,error,*999)

    !Default search method is a linear search
    CALL LIST_SEARCH_LINEAR(A,VALUE,POSITION,err,error,*999)
    
    Call Exits("LIST_SEARCH_C_INT_ARRAY")
    RETURN
999 CALL Errors("LIST_SEARCH_C_INT_ARRAY",err,error)
    Call Exits("LIST_SEARCH_C_INT_ARRAY")
    RETURN 1
  END SUBROUTINE LIST_SEARCH_C_INT_ARRAY
  
  !
  !================================================================================================================================
  !

  !>Searches a single precision real array list A for VALUE. If the search is successful POSITION contains the index of the position of VALUE in the list otherwise POSITION is zero.
  SUBROUTINE LIST_SEARCH_SP_ARRAY(A,VALUE,POSITION,err,error,*)
  
   !Argument variables
    REAL(SP), INTENT(IN) :: A(:) !<The list to search
    REAL(SP), INTENT(IN) :: VALUE !<The value to search for
    INTEGER(INTG), INTENT(OUT) :: POSITION !<On exit, the position of value in the list. If value does not exist in the list the returned position is zero. 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
   
    CALL Enters("LIST_SEARCH_SP_ARRAY",err,error,*999)

    !Default search method is a linear search
    CALL LIST_SEARCH_LINEAR(A,VALUE,POSITION,err,error,*999)    
    
    Call Exits("LIST_SEARCH_SP_ARRAY")
    RETURN
999 CALL Errors("LIST_SEARCH_SP_ARRAY",err,error)
    Call Exits("LIST_SEARCH_SP_ARRAY")
    RETURN 1
  END SUBROUTINE LIST_SEARCH_SP_ARRAY
  
  !
  !================================================================================================================================
  !

  !>Searches a double precision real array list A for VALUE. If the search is successful POSITION contains the index of the position of VALUE in the list otherwise POSITION is zero.
  SUBROUTINE LIST_SEARCH_DP_ARRAY(A,VALUE,POSITION,err,error,*)
  
    !Argument variables
    REAL(DP), INTENT(IN) :: A(:) !<The list to search
    REAL(DP), INTENT(IN) :: VALUE !<The value to search for
    INTEGER(INTG), INTENT(OUT) :: POSITION !<On exit, the position of value in the list. If value does not exist in the list the returned position is zero. 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    
    CALL Enters("LIST_SEARCH_DP_ARRAY",err,error,*999)

    !Default search method is a linear search
    CALL LIST_SEARCH_LINEAR(A,VALUE,POSITION,err,error,*999)    
    
    Call Exits("LIST_SEARCH_DP_ARRAY")
    RETURN
999 CALL Errors("LIST_SEARCH_DP_ARRAY",err,error)
    Call Exits("LIST_SEARCH_DP_ARRAY")
    RETURN 1
  END SUBROUTINE LIST_SEARCH_DP_ARRAY
  
  !
  !================================================================================================================================
  !

  !>Searches an integer array list A for VALUE using the linear search method. If the search is successful POSITION contains the index of the position of VALUE in the list otherwise POSITION is zero.
  SUBROUTINE LIST_SEARCH_LINEAR_INTG_ARRAY(A,VALUE,POSITION,err,error,*)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: A(:) !<The list to search
    INTEGER(INTG), INTENT(IN) :: VALUE !<The value to search for
    INTEGER(INTG), INTENT(OUT) :: POSITION !<On exit, the position of value in the list. If value does not exist in the list the returned position is zero. 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: i
    LOGICAL :: FOUND
    
    CALL Enters("LIST_SEARCH_LINEAR_INTG_ARRAY",err,error,*999)

    FOUND=.FALSE.
    i=1
    DO WHILE(i<=SIZE(A,1).AND..NOT.FOUND)
      IF(A(i)==VALUE) THEN
        FOUND=.TRUE.
      ELSE
        i=i+1
      ENDIF
    ENDDO
    IF(FOUND) THEN
      POSITION=i
    ELSE
      POSITION=0
    ENDIF
    
    Call Exits("LIST_SEARCH_LINEAR_INTG_ARRAY")
    RETURN
999 CALL Errors("LIST_SEARCH_LINEAR_INTG_ARRAY",err,error)
    Call Exits("LIST_SEARCH_LINEAR_INTG_ARRAY")
    RETURN 1
  END SUBROUTINE LIST_SEARCH_LINEAR_INTG_ARRAY
  
  !
  !================================================================================================================================
  !

  !>Searches an integer array list A for VALUE using the linear search method. If the search is successful POSITION contains the index of the position of VALUE in the list otherwise POSITION is zero.
  SUBROUTINE LIST_SEARCH_LINEAR_C_INT_ARRAY(A,VALUE,POSITION,err,error,*)
  
    !Argument variables
    INTEGER(C_INT), INTENT(IN) :: A(:) !<The list to search
    INTEGER(C_INT), INTENT(IN) :: VALUE !<The value to search for
    INTEGER(INTG), INTENT(OUT) :: POSITION !<On exit, the position of value in the list. If value does not exist in the list the returned position is zero. 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: i
    LOGICAL :: FOUND
    
    CALL Enters("LIST_SEARCH_LINEAR_C_INT_ARRAY",err,error,*999)

    FOUND=.FALSE.
    i=1
    DO WHILE(i<=SIZE(A,1).AND..NOT.FOUND)
      IF(A(i)==VALUE) THEN
        FOUND=.TRUE.
      ELSE
        i=i+1
      ENDIF
    ENDDO
    IF(FOUND) THEN
      POSITION=i
    ELSE
      POSITION=0
    ENDIF
    
    Call Exits("LIST_SEARCH_LINEAR_C_INT_ARRAY")
    RETURN
999 CALL Errors("LIST_SEARCH_LINEAR_C_INT_ARRAY",err,error)
    Call Exits("LIST_SEARCH_LINEAR_C_INT_ARRAY")
    RETURN 1
  END SUBROUTINE LIST_SEARCH_LINEAR_C_INT_ARRAY
  
  !
  !================================================================================================================================
  !

  !>Searches a single precision real array list A for VALUE using the linear search method. If the search is successful POSITION contains the index of the position of VALUE in the list otherwise POSITION is zero.
  SUBROUTINE LIST_SEARCH_LINEAR_SP_ARRAY(A,VALUE,POSITION,err,error,*)
  
    !Argument variables
    REAL(SP), INTENT(IN) :: A(:) !<The list to search
    REAL(SP), INTENT(IN) :: VALUE !<The value to search for
    INTEGER(INTG), INTENT(OUT) :: POSITION !<On exit, the position of value in the list. If value does not exist in the list the returned position is zero. 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: i
    LOGICAL :: FOUND
    
    CALL Enters("LIST_SEARCH_LINEAR_SP_ARRAY",err,error,*999)

    FOUND=.FALSE.
    i=1
    DO WHILE(i<=SIZE(A,1).AND..NOT.FOUND)
      IF(ABS(A(i)-VALUE)<ZERO_TOLERANCE_SP) THEN
        FOUND=.TRUE.
      ELSE
        i=i+1
      ENDIF
    ENDDO
    IF(FOUND) THEN
      POSITION=i
    ELSE
      POSITION=0
    ENDIF
    
    Call Exits("LIST_SEARCH_LINEAR_SP_ARRAY")
    RETURN
999 CALL Errors("LIST_SEARCH_LINEAR_SP_ARRAY",err,error)
    Call Exits("LIST_SEARCH_LINEAR_SP_ARRAY")
    RETURN 1
  END SUBROUTINE LIST_SEARCH_LINEAR_SP_ARRAY
  
  !
  !================================================================================================================================
  !

  !>Searches a double precision real array list A for VALUE using the linear search method. If the search is successful POSITION contains the index of the position of VALUE in the list otherwise POSITION is zero.
  SUBROUTINE LIST_SEARCH_LINEAR_DP_ARRAY(A,VALUE,POSITION,err,error,*)
  
    !Argument variables
    REAL(DP), INTENT(IN) :: A(:) !<The list to search
    REAL(DP), INTENT(IN) :: VALUE !<The value to search for
    INTEGER(INTG), INTENT(OUT) :: POSITION !<On exit, the position of value in the list. If value does not exist in the list the returned position is zero. 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: i
    LOGICAL :: FOUND
    
    CALL Enters("LIST_SEARCH_LINEAR_DP_ARRAY",err,error,*999)

    FOUND=.FALSE.
    i=1
    DO WHILE(i<=SIZE(A,1).AND..NOT.FOUND)
      IF(ABS(A(i)-VALUE)<ZERO_TOLERANCE) THEN
        FOUND=.TRUE.
      ELSE
        i=i+1
      ENDIF
    ENDDO
    IF(FOUND) THEN
      POSITION=i
    ELSE
      POSITION=0
    ENDIF
    
    Call Exits("LIST_SEARCH_LINEAR_DP_ARRAY")
    RETURN
999 CALL Errors("LIST_SEARCH_LINEAR_DP_ARRAY",err,error)
    Call Exits("LIST_SEARCH_LINEAR_DP_ARRAY")
    RETURN 1
  END SUBROUTINE LIST_SEARCH_LINEAR_DP_ARRAY
  
  !
  !================================================================================================================================
  !

  !>Sorts a list of into ascending order.
  SUBROUTINE List_SortList(list,err,error,*)
  
    !Argument variables
    TYPE(LIST_TYPE), POINTER, INTENT(INOUT) :: LIST !<The list to sort
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localError
    
    CALL Enters("List_SortList",err,error,*999)

    IF(ASSOCIATED(list)) THEN
      IF(list%LIST_FINISHED) THEN
        SELECT CASE(list%SORT_METHOD)
        CASE(LIST_BUBBLE_SORT_METHOD)
          IF(list%DATA_DIMENSION==1) THEN
            SELECT CASE(list%DATA_TYPE)
            CASE(LIST_INTG_TYPE)
              CALL LIST_SORT_BUBBLE_INTG1_ARRAY(list%LIST_INTG(1:list%NUMBER_IN_LIST),ERR,ERROR,*999)
            CASE(LIST_SP_TYPE)
              CALL LIST_SORT_BUBBLE_SP1_ARRAY(list%LIST_SP(1:list%NUMBER_IN_LIST),err,error,*999)              
            CASE(LIST_DP_TYPE)
              CALL LIST_SORT_BUBBLE_DP1_ARRAY(list%LIST_DP(1:list%NUMBER_IN_LIST),err,error,*999) 
            CASE DEFAULT
              localError="The list data type of "//TRIM(NumberToVString(list%DATA_TYPE,"*",err,error))//" is invalid."
              CALL FlagError(localError,err,error,*999)
            END SELECT            
          ELSE
            SELECT CASE(list%DATA_TYPE)
            CASE(LIST_INTG_TYPE)
              CALL LIST_SORT_BUBBLE_INTG2_ARRAY(list%LIST_INTG2(:,1:list%NUMBER_IN_LIST),list%KEY_DIMENSION, &
                & err,error,*999)
            CASE(LIST_SP_TYPE)
              CALL LIST_SORT_BUBBLE_SP2_ARRAY(list%LIST_SP2(:,1:list%NUMBER_IN_LIST),list%KEY_DIMENSION, &
                & err,error,*999)              
            CASE(LIST_DP_TYPE)
              CALL LIST_SORT_BUBBLE_DP2_ARRAY(list%LIST_DP2(:,1:list%NUMBER_IN_LIST),list%KEY_DIMENSION, &
                & err,error,*999)                            
            CASE DEFAULT
              localError="The list data type of "//TRIM(NumberToVString(list%DATA_TYPE,"*",err,error))//" is invalid."
              CALL FlagError(localError,err,error,*999)
            END SELECT            
          ENDIF
        CASE(LIST_SHELL_SORT_METHOD)
          IF(list%DATA_DIMENSION==1) THEN
            SELECT CASE(list%DATA_TYPE)
            CASE(LIST_INTG_TYPE)
              CALL LIST_SORT_SHELL_INTG1_ARRAY(list%LIST_INTG(1:list%NUMBER_IN_LIST),err,error,*999)
           CASE(LIST_SP_TYPE)
              CALL LIST_SORT_SHELL_SP1_ARRAY(list%LIST_SP(1:list%NUMBER_IN_LIST),err,error,*999)              
            CASE(LIST_DP_TYPE)
              CALL LIST_SORT_SHELL_DP1_ARRAY(list%LIST_DP(1:list%NUMBER_IN_LIST),err,error,*999)
            CASE DEFAULT
              localError="The list data type of "//TRIM(NumberToVString(list%DATA_TYPE,"*",err,error))//" is invalid."
              CALL FlagError(localError,err,error,*999)
            END SELECT            
          ELSE
            SELECT CASE(list%DATA_TYPE)
            CASE(LIST_INTG_TYPE)
              CALL LIST_SORT_SHELL_INTG2_ARRAY(list%LIST_INTG2(:,1:list%NUMBER_IN_LIST),list%KEY_DIMENSION, &
                & err,error,*999)
            CASE(LIST_SP_TYPE)
              CALL LIST_SORT_SHELL_SP2_ARRAY(list%LIST_SP2(:,1:list%NUMBER_IN_LIST),list%KEY_DIMENSION, &
                & err,error,*999)              
            CASE(LIST_DP_TYPE)
              CALL LIST_SORT_SHELL_DP2_ARRAY(list%LIST_DP2(:,1:list%NUMBER_IN_LIST),list%KEY_DIMENSION, &
                & err,error,*999)                            
            CASE DEFAULT
              localError="The list data type of "//TRIM(NumberToVString(list%DATA_TYPE,"*",err,error))//" is invalid."
              CALL FlagError(localError,err,error,*999)
            END SELECT            
          ENDIF
        CASE(LIST_HEAP_SORT_METHOD)
          IF(list%DATA_DIMENSION==1) THEN
            SELECT CASE(list%DATA_TYPE)
            CASE(LIST_INTG_TYPE)
              CALL LIST_SORT_HEAP_INTG1_ARRAY(list%LIST_INTG(1:list%NUMBER_IN_LIST),err,error,*999)
            CASE(LIST_SP_TYPE)
              CALL LIST_SORT_HEAP_SP1_ARRAY(list%LIST_SP(1:list%NUMBER_IN_LIST),err,error,*999)              
            CASE(LIST_DP_TYPE)
              CALL LIST_SORT_HEAP_DP1_ARRAY(list%LIST_DP(1:list%NUMBER_IN_LIST),err,error,*999)                            
            CASE DEFAULT
              localError="The list data type of "//TRIM(NumberToVString(list%DATA_TYPE,"*",err,error))//" is invalid."
              CALL FlagError(localError,err,error,*999)
            END SELECT            
          ELSE
            SELECT CASE(list%DATA_TYPE)
            CASE(LIST_INTG_TYPE)
              CALL LIST_SORT_HEAP_INTG2_ARRAY(list%LIST_INTG2(:,1:list%NUMBER_IN_LIST),list%KEY_DIMENSION, &
                & err,error,*999)
            CASE(LIST_SP_TYPE)
              CALL LIST_SORT_HEAP_SP2_ARRAY(list%LIST_SP2(:,1:list%NUMBER_IN_LIST),list%KEY_DIMENSION, &
                & err,error,*999)              
            CASE(LIST_DP_TYPE)
              CALL LIST_SORT_HEAP_DP2_ARRAY(list%LIST_DP2(:,1:list%NUMBER_IN_LIST),list%KEY_DIMENSION, &
                & err,error,*999)                            
            CASE DEFAULT
              localError="The list data type of "//TRIM(NumberToVString(list%DATA_TYPE,"*",err,error))//" is invalid."
              CALL FlagError(localError,err,error,*999)
            END SELECT            
          ENDIF
        CASE DEFAULT
          localError="The list sort method of "//TRIM(NumberToVString(list%SORT_METHOD,"*",err,error))//" is invlaid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ELSE
        CALL FlagError("List has not been finished.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("List is not associated.",err,error,*999)
    ENDIF

    Call Exits("List_SortList")
    RETURN
999 CALL Errors("List_SortList",err,error)
    Call Exits("List_SortList")
    RETURN 1
  END SUBROUTINE List_SortList
  
  !
  !================================================================================================================================
  !

  !>Sorts a list of into ascending order.
  SUBROUTINE LIST_SORT_LIST(LIST,err,error,*)
  
    !Argument variables
    TYPE(LIST_TYPE), POINTER, INTENT(INOUT) :: LIST !<The list to sort
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables

    CALL List_SortList(LIST,err,error,*999)
    
    RETURN
999 RETURN 1
  END SUBROUTINE LIST_SORT_LIST
  
  !
  !================================================================================================================================
  !

  !>Sorts an integer array list of data dimension 1 into ascending order.
  SUBROUTINE LIST_SORT_INTG1_ARRAY(A,err,error,*)
  
    !Argument variables
    INTEGER(INTG), INTENT(INOUT) :: A(:) !<The list to sort
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    
    CALL Enters("LIST_SORT_INTG1_ARRAY",err,error,*999)

    !Default sort method is a heap sort
    CALL LIST_SORT_HEAP(A,err,error,*999)    

    Call Exits("LIST_SORT_INTG1_ARRAY")
    RETURN
999 CALL Errors("LIST_SORT_INTG1_ARRAY",err,error)
    Call Exits("LIST_SORT_INTG1_ARRAY")
    RETURN 1
  END SUBROUTINE LIST_SORT_INTG1_ARRAY
  
   !
  !================================================================================================================================
  !

  !>Sorts an integer array list of data dimension > 1 into ascending order.
  SUBROUTINE LIST_SORT_INTG2_ARRAY(A,KEY_DIMENSION,err,error,*)
  
    !Argument variables
    INTEGER(INTG), INTENT(INOUT) :: A(:,:) !<The list to sort
    INTEGER(INTG), INTENT(IN) :: KEY_DIMENSION !<The key dimension to sort on
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    
    CALL Enters("LIST_SORT_INTG2_ARRAY",err,error,*999)

    !Default sort method is a heap sort
    CALL LIST_SORT_HEAP(A,KEY_DIMENSION,err,error,*999)    

    Call Exits("LIST_SORT_INTG2_ARRAY")
    RETURN
999 CALL Errors("LIST_SORT_INTG2_ARRAY",err,error)
    Call Exits("LIST_SORT_INTG2_ARRAY")
    RETURN 1
  END SUBROUTINE LIST_SORT_INTG2_ARRAY
  
  !
  !================================================================================================================================
  !

  !>Sorts an integer array list of data dimension 1 into ascending order.
  SUBROUTINE LIST_SORT_C_INT1_ARRAY(A,err,error,*)
  
    !Argument variables
    INTEGER(C_INT), INTENT(INOUT) :: A(:) !<The list to sort
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    
    CALL Enters("LIST_SORT_C_INT1_ARRAY",err,error,*999)

    !Default sort method is a heap sort
    CALL LIST_SORT_HEAP(A,err,error,*999)    

    Call Exits("LIST_SORT_C_INT1_ARRAY")
    RETURN
999 CALL Errors("LIST_SORT_C_INT1_ARRAY",err,error)
    Call Exits("LIST_SORT_C_INT1_ARRAY")
    RETURN 1
  END SUBROUTINE LIST_SORT_C_INT1_ARRAY
  
   !
  !================================================================================================================================
  !

  !>Sorts an integer array list of data dimension > 1 into ascending order.
  SUBROUTINE LIST_SORT_C_INT2_ARRAY(A,KEY_DIMENSION,err,error,*)
  
    !Argument variables
    INTEGER(C_INT), INTENT(INOUT) :: A(:,:) !<The list to sort
    INTEGER(INTG), INTENT(IN) :: KEY_DIMENSION !<The key dimension to sort on
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    
    CALL Enters("LIST_SORT_C_INT2_ARRAY",err,error,*999)

    !Default sort method is a heap sort
    CALL LIST_SORT_HEAP(A,KEY_DIMENSION,err,error,*999)    

    Call Exits("LIST_SORT_C_INT2_ARRAY")
    RETURN
999 CALL Errors("LIST_SORT_C_INT2_ARRAY",err,error)
    Call Exits("LIST_SORT_C_INT2_ARRAY")
    RETURN 1
  END SUBROUTINE LIST_SORT_C_INT2_ARRAY
  
  !
  !================================================================================================================================
  !

  !>Sorts an single precision array list of data dimension 1 into ascending order.
  SUBROUTINE LIST_SORT_SP1_ARRAY(A,err,error,*)
      
    !Argument variables
    REAL(SP), INTENT(INOUT) :: A(:) !<The list to sort
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    
    CALL Enters("LIST_SORT_SP1_ARRAY",err,error,*999)

    !Default sort method is a heap sort
    CALL LIST_SORT_HEAP(A,err,error,*999)    

    Call Exits("LIST_SORT_SP1_ARRAY")
    RETURN
999 CALL Errors("LIST_SORT_SP1_ARRAY",err,error)
    Call Exits("LIST_SORT_SP1_ARRAY")
    RETURN 1
  END SUBROUTINE LIST_SORT_SP1_ARRAY
  
  !
  !================================================================================================================================
  !

  !>Sorts an single precision array list of data dimension > 1 into ascending order.
  SUBROUTINE LIST_SORT_SP2_ARRAY(A,KEY_DIMENSION,err,error,*)
      
    !Argument variables
    REAL(SP), INTENT(INOUT) :: A(:,:) !<The list to sort
    INTEGER(INTG), INTENT(IN) :: KEY_DIMENSION !<The key dimension to sort the list on.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    
    CALL Enters("LIST_SORT_SP2_ARRAY",err,error,*999)

    !Default sort method is a heap sort
    CALL LIST_SORT_HEAP(A,KEY_DIMENSION,err,error,*999)    

    Call Exits("LIST_SORT_SP2_ARRAY")
    RETURN
999 CALL Errors("LIST_SORT_SP2_ARRAY",err,error)
    Call Exits("LIST_SORT_SP2_ARRAY")
    RETURN 1
  END SUBROUTINE LIST_SORT_SP2_ARRAY
  
  !
  !================================================================================================================================
  !

  !>Sorts an double precision array list of data dimension 1 into ascending order.
  SUBROUTINE LIST_SORT_DP1_ARRAY(A,err,error,*)
  
    !Argument variables
    REAL(DP), INTENT(INOUT) :: A(:) !<The list to sort
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
     
    CALL Enters("LIST_SORT_DP1_ARRAY",err,error,*999)

    !Default sort method is a heap sort
    CALL LIST_SORT_HEAP(A,err,error,*999)    

    Call Exits("LIST_SORT_DP1_ARRAY")
    RETURN
999 CALL Errors("LIST_SORT_DP1_ARRAY",err,error)
    Call Exits("LIST_SORT_DP1_ARRAY")
    RETURN 1
  END SUBROUTINE LIST_SORT_DP1_ARRAY
  
  !
  !================================================================================================================================
  !

  !>Sorts an double precision array list of data dimension > 1 into ascending order.
  SUBROUTINE LIST_SORT_DP2_ARRAY(A,KEY_DIMENSION,err,error,*)
  
    !Argument variables
    REAL(DP), INTENT(INOUT) :: A(:,:) !<The list to sort
    INTEGER(INTG), INTENT(IN) :: KEY_DIMENSION !<The key dimension to sort on.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
     
    CALL Enters("LIST_SORT_DP2_ARRAY",err,error,*999)

    !Default sort method is a heap sort
    CALL LIST_SORT_HEAP(A,KEY_DIMENSION,err,error,*999)    

    Call Exits("LIST_SORT_DP2_ARRAY")
    RETURN
999 CALL Errors("LIST_SORT_DP2_ARRAY",err,error)
    Call Exits("LIST_SORT_DP2_ARRAY")
    RETURN 1
  END SUBROUTINE LIST_SORT_DP2_ARRAY
  
  !
  !================================================================================================================================
  !

  !>BUBBLE_SORT_INTG performs a bubble sort on an integer array of data dimension 1 list
  SUBROUTINE LIST_SORT_BUBBLE_INTG1_ARRAY(A,err,error,*)
  
    !Argument variables
    INTEGER(INTG), INTENT(INOUT) :: A(:) !<The list to sort
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: FLAG,i,j,k,VALUE
    
    CALL Enters("LIST_SORT_BUBBLE_INTG1_ARRAY",err,error,*999)

    IF(SIZE(A,1)>1) THEN
      FLAG=SIZE(A,1)
      DO i=1,SIZE(A,1)
        k=FLAG-1
        FLAG=0
        DO j=1,k
          IF(A(j)>A(j+1)) THEN
            VALUE=A(j)
            A(j)=A(j+1)
            A(j+1)=VALUE
            FLAG=j
          ENDIF
        ENDDO
        IF(FLAG==0) EXIT
      ENDDO
    ENDIF

    Call Exits("LIST_SORT_BUBBLE_INTG1_ARRAY")
    RETURN
999 CALL Errors("LIST_SORT_BUBBLE_INTG1_ARRAY",err,error)
    Call Exits("LIST_SORT_BUBBLE_INTG1_ARRAY")
    RETURN 1
  END SUBROUTINE LIST_SORT_BUBBLE_INTG1_ARRAY
  
  !
  !================================================================================================================================
  !

  !>BUBBLE_SORT_INTG performs a bubble sort on an integer array of data dimension > 1 list
  SUBROUTINE LIST_SORT_BUBBLE_INTG2_ARRAY(A,KEY_DIMENSION,err,error,*)
  
    !Argument variables
    INTEGER(INTG), INTENT(INOUT) :: A(:,:) !<The list to sort
    INTEGER(INTG), INTENT(IN) :: KEY_DIMENSION !<The key dimension of A to do the sort on
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: FLAG,i,j,k,VALUE(SIZE(A,1))
    TYPE(VARYING_STRING) :: localError
    
    CALL Enters("LIST_SORT_BUBBLE_INTG2_ARRAY",err,error,*999)

    IF(KEY_DIMENSION>0.AND.KEY_DIMENSION<=SIZE(A,1)) THEN
      IF(SIZE(A,2)>1) THEN
        FLAG=SIZE(A,2)
        DO i=1,SIZE(A,2)
          k=FLAG-1
          FLAG=0
          DO j=1,k
            IF(A(KEY_DIMENSION,j)>A(KEY_DIMENSION,j+1)) THEN
              VALUE=A(:,j)
              A(:,j)=A(:,j+1)
              A(:,j+1)=VALUE
              FLAG=j
            ENDIF
          ENDDO
          IF(FLAG==0) EXIT
        ENDDO
      ENDIF
    ELSE
      localError="The specified key dimension of "//TRIM(NumberToVString(KEY_DIMENSION,"*",err,error))// &
        & " is invalid. The key dimension must be > 0 and <= "//TRIM(NumberToVString(SIZE(A,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    Call Exits("LIST_SORT_BUBBLE_INTG2_ARRAY")
    RETURN
999 CALL Errors("LIST_SORT_BUBBLE_INTG2_ARRAY",err,error)
    Call Exits("LIST_SORT_BUBBLE_INTG2_ARRAY")
    RETURN 1
  END SUBROUTINE LIST_SORT_BUBBLE_INTG2_ARRAY
  
  !
  !================================================================================================================================
  !

  !>BUBBLE_SORT_C_INT performs a bubble sort on an integer array of data dimension 1 list
  SUBROUTINE LIST_SORT_BUBBLE_C_INT1_ARRAY(A,err,error,*)
  
    !Argument variables
    INTEGER(C_INT), INTENT(INOUT) :: A(:) !<The list to sort
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: FLAG,i,j,k
    INTEGER(C_INT) :: VALUE
    
    CALL Enters("LIST_SORT_BUBBLE_C_INT1_ARRAY",err,error,*999)

    IF(SIZE(A,1)>1) THEN
      FLAG=SIZE(A,1)
      DO i=1,SIZE(A,1)
        k=FLAG-1
        FLAG=0
        DO j=1,k
          IF(A(j)>A(j+1)) THEN
            VALUE=A(j)
            A(j)=A(j+1)
            A(j+1)=VALUE
            FLAG=j
          ENDIF
        ENDDO
        IF(FLAG==0) EXIT
      ENDDO
    ENDIF

    Call Exits("LIST_SORT_BUBBLE_C_INT1_ARRAY")
    RETURN
999 CALL Errors("LIST_SORT_BUBBLE_C_INT1_ARRAY",err,error)
    Call Exits("LIST_SORT_BUBBLE_C_INT1_ARRAY")
    RETURN 1
  END SUBROUTINE LIST_SORT_BUBBLE_C_INT1_ARRAY
  
  !
  !================================================================================================================================
  !

  !>BUBBLE_SORT_C_INT performs a bubble sort on an integer array of data dimension > 1 list
  SUBROUTINE LIST_SORT_BUBBLE_C_INT2_ARRAY(A,KEY_DIMENSION,err,error,*)
  
    !Argument variables
    INTEGER(C_INT), INTENT(INOUT) :: A(:,:) !<The list to sort
    INTEGER(INTG), INTENT(IN) :: KEY_DIMENSION !<The key dimension of A to do the sort on
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: FLAG,i,j,k
    INTEGER(C_INT) :: VALUE(SIZE(A,1))
    TYPE(VARYING_STRING) :: localError
    
    CALL Enters("LIST_SORT_BUBBLE_C_INT2_ARRAY",err,error,*999)

    IF(KEY_DIMENSION>0.AND.KEY_DIMENSION<=SIZE(A,1)) THEN
      IF(SIZE(A,2)>1) THEN
        FLAG=SIZE(A,2)
        DO i=1,SIZE(A,2)
          k=FLAG-1
          FLAG=0
          DO j=1,k
            IF(A(KEY_DIMENSION,j)>A(KEY_DIMENSION,j+1)) THEN
              VALUE=A(:,j)
              A(:,j)=A(:,j+1)
              A(:,j+1)=VALUE
              FLAG=j
            ENDIF
          ENDDO
          IF(FLAG==0) EXIT
        ENDDO
      ENDIF
    ELSE
      localError="The specified key dimension of "//TRIM(NumberToVString(KEY_DIMENSION,"*",err,error))// &
        & " is invalid. The key dimension must be > 0 and <= "//TRIM(NumberToVString(SIZE(A,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    Call Exits("LIST_SORT_BUBBLE_C_INT2_ARRAY")
    RETURN
999 CALL Errors("LIST_SORT_BUBBLE_C_INT2_ARRAY",err,error)
    Call Exits("LIST_SORT_BUBBLE_C_INT2_ARRAY")
    RETURN 1
  END SUBROUTINE LIST_SORT_BUBBLE_C_INT2_ARRAY
  
  !
  !================================================================================================================================
  !

  !>BUBBLE_SORT_SP performs a bubble sort on a single precision array of data dimension 1 list
  SUBROUTINE LIST_SORT_BUBBLE_SP1_ARRAY(A,err,error,*)
  
    !Argument variables
    REAL(SP), INTENT(INOUT) :: A(:) !<The list to sort
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: FLAG,i,j,k
    REAL(SP) :: VALUE
    
    CALL Enters("LIST_SORT_BUBBLE_SP1_ARRAY",err,error,*999)

    IF(SIZE(A,1)>1) THEN
      FLAG=SIZE(A,1)
      DO i=1,SIZE(A,1)
        k=FLAG-1
        FLAG=0
        DO j=1,k
          IF(A(j)>A(j+1)) THEN
            VALUE=A(j)
            A(j)=A(j+1)
            A(j+1)=VALUE
            FLAG=j
          ENDIF
        ENDDO
        IF(FLAG==0) EXIT
      ENDDO
    ENDIF

    Call Exits("LIST_SORT_BUBBLE_SP1_ARRAY")
    RETURN
999 CALL Errors("LIST_SORT_BUBBLE_SP1_ARRAY",err,error)
    Call Exits("LIST_SORT_BUBBLE_SP1_ARRAY")
    RETURN 1
  END SUBROUTINE LIST_SORT_BUBBLE_SP1_ARRAY
  
   !
  !================================================================================================================================
  !

  !>BUBBLE_SORT_SP performs a bubble sort on a single precision array of data dimension > 1 list
  SUBROUTINE LIST_SORT_BUBBLE_SP2_ARRAY(A,KEY_DIMENSION,err,error,*)
  
    !Argument variables
    REAL(SP), INTENT(INOUT) :: A(:,:) !<The list to sort
    INTEGER(INTG), INTENT(IN) :: KEY_DIMENSION !<The key dimension of A to do the sort on
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: FLAG,i,j,k
    REAL(SP) :: VALUE(SIZE(A,1))
    TYPE(VARYING_STRING) :: localError
     
    CALL Enters("LIST_SORT_BUBBLE_SP2_ARRAY",err,error,*999)

    IF(KEY_DIMENSION>0.AND.KEY_DIMENSION<=SIZE(A,1)) THEN
      IF(SIZE(A,2)>1) THEN
        FLAG=SIZE(A,2)
        DO i=1,SIZE(A,2)
          k=FLAG-1
          FLAG=0
          DO j=1,k
            IF(A(KEY_DIMENSION,j)>A(KEY_DIMENSION,j+1)) THEN
              VALUE=A(:,j)
              A(:,j)=A(:,j+1)
              A(:,j+1)=VALUE
              FLAG=j
            ENDIF
          ENDDO
          IF(FLAG==0) EXIT
        ENDDO
      ENDIF
    ELSE
      localError="The specified key dimension of "//TRIM(NumberToVString(KEY_DIMENSION,"*",err,error))// &
        & " is invalid. The key dimension must be > 0 and <= "//TRIM(NumberToVString(SIZE(A,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    Call Exits("LIST_SORT_BUBBLE_SP2_ARRAY")
    RETURN
999 CALL Errors("LIST_SORT_BUBBLE_SP2_ARRAY",err,error)
    Call Exits("LIST_SORT_BUBBLE_SP2_ARRAY")
    RETURN 1
  END SUBROUTINE LIST_SORT_BUBBLE_SP2_ARRAY
  
  !
  !================================================================================================================================
  !

  !>BUBBLE_SORT_DP performs a bubble sort on a double precision of data dimension 1 list
  SUBROUTINE LIST_SORT_BUBBLE_DP1_ARRAY(A,err,error,*)
  
    !Argument variables
    REAL(DP), INTENT(INOUT) :: A(:) !<The list to sort
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: FLAG,i,j,k
    REAL(DP) :: VALUE
    
    CALL Enters("LIST_SORT_BUBBLE_DP1_ARRAY",err,error,*999)

    IF(SIZE(A,1)>1) THEN
      FLAG=SIZE(A,1)
      DO i=1,SIZE(A,1)
        k=FLAG-1
        FLAG=0
        DO j=1,k
          IF(A(j)>A(j+1)) THEN
            VALUE=A(j)
            A(j)=A(j+1)
            A(j+1)=VALUE
            FLAG=j
          ENDIF
        ENDDO
        IF(FLAG==0) EXIT
      ENDDO
    ENDIF

    Call Exits("LIST_SORT_BUBBLE_DP1_ARRAY")
    RETURN
999 CALL Errors("LIST_SORT_BUBBLE_DP1_ARRAY",err,error)
    Call Exits("LIST_SORT_BUBBLE_DP1_ARRAY")
    RETURN 1
  END SUBROUTINE LIST_SORT_BUBBLE_DP1_ARRAY
  
  !
  !================================================================================================================================
  !

  !>BUBBLE_SORT_DP performs a bubble sort on a double precision of data dimension > 1 list
  SUBROUTINE LIST_SORT_BUBBLE_DP2_ARRAY(A,KEY_DIMENSION,err,error,*)
  
    !Argument variables
    REAL(DP), INTENT(INOUT) :: A(:,:) !<The list to sort
    INTEGER(INTG), INTENT(IN) :: KEY_DIMENSION !<The key dimension of A to do the sort on
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: FLAG,i,j,k
    REAL(DP) :: VALUE(SIZE(A,1))
    TYPE(VARYING_STRING) :: localError
    
    CALL Enters("LIST_SORT_BUBBLE_DP2_ARRAY",err,error,*999)

    IF(KEY_DIMENSION>0.AND.KEY_DIMENSION<=SIZE(A,1)) THEN
      IF(SIZE(A,2)>1) THEN
        FLAG=SIZE(A,2)
        DO i=1,SIZE(A,2)
          k=FLAG-1
          FLAG=0
          DO j=1,k
            IF(A(KEY_DIMENSION,j)>A(KEY_DIMENSION,j+1)) THEN
              VALUE=A(:,j)
              A(:,j)=A(:,j+1)
              A(:,j+1)=VALUE
              FLAG=j
            ENDIF
          ENDDO
          IF(FLAG==0) EXIT
        ENDDO
      ENDIF
    ELSE
      localError="The specified key dimension of "//TRIM(NumberToVString(KEY_DIMENSION,"*",err,error))// &
        & " is invalid. The key dimension must be > 0 and <= "//TRIM(NumberToVString(SIZE(A,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    Call Exits("LIST_SORT_BUBBLE_DP2_ARRAY")
    RETURN
999 CALL Errors("LIST_SORT_BUBBLE_DP2_ARRAY",err,error)
    Call Exits("LIST_SORT_BUBBLE_DP2_ARRAY")
    RETURN 1
  END SUBROUTINE LIST_SORT_BUBBLE_DP2_ARRAY
  
  !
  !================================================================================================================================
  !

  !>Sorts an integer array of data dimension 1 list into assending order using the heap sort method.
  SUBROUTINE LIST_SORT_HEAP_INTG1_ARRAY(A,err,error,*)
  
    !Argument variables
    INTEGER(INTG), INTENT(INOUT) :: A(:) !<The list to sort
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: I,IVALUE,J,L,VALUE
    
    CALL Enters("LIST_SORT_HEAP_INTG1_ARRAY",err,error,*999)

    IF(SIZE(A,1)>1) THEN      
      L=SIZE(A,1)/2+1
      IVALUE=SIZE(A,1)
      DO 
        IF(L>1) THEN
          L=L-1
          VALUE=A(L)
        ELSE
          VALUE=A(IVALUE)
          A(IVALUE)=A(1)
          IVALUE=IVALUE-1
          IF(IVALUE==1) THEN
            A(1)=VALUE
            EXIT
          ENDIF
        ENDIF
        I=L
        J=L+L
        DO WHILE(J<=IVALUE)
          IF(J<IVALUE) THEN
            IF(A(J)<A(J+1)) J=J+1
          ENDIF
          IF(VALUE<A(J)) THEN
            A(I)=A(J)
            I=J
            J=J+J
          ELSE
            J=IVALUE+1
          ENDIF
        ENDDO
        A(I)=VALUE
      ENDDO
    ENDIF

    Call Exits("LIST_SORT_HEAP_INTG1_ARRAY")
    RETURN
999 CALL Errors("LIST_SORT_HEAP_INTG1_ARRAY",err,error)
    Call Exits("LIST_SORT_HEAP_INTG1_ARRAY")
    RETURN 1
  END SUBROUTINE LIST_SORT_HEAP_INTG1_ARRAY
  
  !
  !================================================================================================================================
  !

  !>Sorts an integer array of data dimension > 1 list into assending order using the heap sort method.
  SUBROUTINE LIST_SORT_HEAP_INTG2_ARRAY(A,KEY_DIMENSION,err,error,*)
  
    !Argument variables
    INTEGER(INTG), INTENT(INOUT) :: A(:,:) !<The list to sort
    INTEGER(INTG), INTENT(IN) :: KEY_DIMENSION !<The key dimension of A to do the sort on
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: I,IVALUE,J,L,VALUE(SIZE(A,1))
    TYPE(VARYING_STRING) :: localError
    
    CALL Enters("LIST_SORT_HEAP_INTG2_ARRAY",err,error,*999)
    
    IF(KEY_DIMENSION>0.AND.KEY_DIMENSION<=SIZE(A,1)) THEN
      IF(SIZE(A,2)>1) THEN      
        L=SIZE(A,2)/2+1
        IVALUE=SIZE(A,2)
        DO 
          IF(L>1) THEN
            L=L-1
            VALUE=A(:,L)
          ELSE
            VALUE=A(:,IVALUE)
            A(:,IVALUE)=A(:,1)
            IVALUE=IVALUE-1
            IF(IVALUE==1) THEN
              A(:,1)=VALUE
              EXIT
            ENDIF
          ENDIF
          I=L
          J=L+L
          DO WHILE(J<=IVALUE)
            IF(J<IVALUE) THEN
              IF(A(KEY_DIMENSION,J)<A(KEY_DIMENSION,J+1)) J=J+1
            ENDIF
            IF(VALUE(KEY_DIMENSION)<A(KEY_DIMENSION,J)) THEN
              A(:,I)=A(:,J)
              I=J
              J=J+J
            ELSE
              J=IVALUE+1
            ENDIF
          ENDDO
          A(:,I)=VALUE
        ENDDO
      ENDIF
    ELSE
      localError="The specified key dimension of "//TRIM(NumberToVString(KEY_DIMENSION,"*",err,error))// &
        & " is invalid. The key dimension must be > 0 and <= "//TRIM(NumberToVString(SIZE(A,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    Call Exits("LIST_SORT_HEAP_INTG2_ARRAY")
    RETURN
999 CALL Errors("LIST_SORT_HEAP_INTG2_ARRAY",err,error)
    Call Exits("LIST_SORT_HEAP_INTG2_ARRAY")
    RETURN 1
  END SUBROUTINE LIST_SORT_HEAP_INTG2_ARRAY
  
  !
  !================================================================================================================================
  !

  !>Sorts an integer array of data dimension 1 list into assending order using the heap sort method.
  SUBROUTINE LIST_SORT_HEAP_C_INT1_ARRAY(A,err,error,*)
  
    !Argument variables
    INTEGER(C_INT), INTENT(INOUT) :: A(:) !<The list to sort
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: I,J,L
    INTEGER(C_INT) :: IVALUE,VALUE
    
    CALL Enters("LIST_SORT_HEAP_C_INT1_ARRAY",err,error,*999)

    IF(SIZE(A,1)>1) THEN      
      L=SIZE(A,1)/2+1
      IVALUE=SIZE(A,1)
      DO 
        IF(L>1) THEN
          L=L-1
          VALUE=A(L)
        ELSE
          VALUE=A(IVALUE)
          A(IVALUE)=A(1)
          IVALUE=IVALUE-1
          IF(IVALUE==1) THEN
            A(1)=VALUE
            EXIT
          ENDIF
        ENDIF
        I=L
        J=L+L
        DO WHILE(J<=IVALUE)
          IF(J<IVALUE) THEN
            IF(A(J)<A(J+1)) J=J+1
          ENDIF
          IF(VALUE<A(J)) THEN
            A(I)=A(J)
            I=J
            J=J+J
          ELSE
            J=IVALUE+1
          ENDIF
        ENDDO
        A(I)=VALUE
      ENDDO
    ENDIF

    Call Exits("LIST_SORT_HEAP_C_INT1_ARRAY")
    RETURN
999 CALL Errors("LIST_SORT_HEAP_C_INT1_ARRAY",err,error)
    Call Exits("LIST_SORT_HEAP_C_INT1_ARRAY")
    RETURN 1
  END SUBROUTINE LIST_SORT_HEAP_C_INT1_ARRAY
  
  !
  !================================================================================================================================
  !

  !>Sorts an integer array of data dimension > 1 list into assending order using the heap sort method.
  SUBROUTINE LIST_SORT_HEAP_C_INT2_ARRAY(A,KEY_DIMENSION,err,error,*)
  
    !Argument variables
    INTEGER(C_INT), INTENT(INOUT) :: A(:,:) !<The list to sort
    INTEGER(INTG), INTENT(IN) :: KEY_DIMENSION !<The key dimension of A to do the sort on
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: I,J,L
    INTEGER(C_INT) :: IVALUE,VALUE(SIZE(A,1))
    TYPE(VARYING_STRING) :: localError
    
    CALL Enters("LIST_SORT_HEAP_C_INT2_ARRAY",err,error,*999)
    
    IF(KEY_DIMENSION>0.AND.KEY_DIMENSION<=SIZE(A,1)) THEN
      IF(SIZE(A,2)>1) THEN      
        L=SIZE(A,2)/2+1
        IVALUE=SIZE(A,2)
        DO 
          IF(L>1) THEN
            L=L-1
            VALUE=A(:,L)
          ELSE
            VALUE=A(:,IVALUE)
            A(:,IVALUE)=A(:,1)
            IVALUE=IVALUE-1
            IF(IVALUE==1) THEN
              A(:,1)=VALUE
              EXIT
            ENDIF
          ENDIF
          I=L
          J=L+L
          DO WHILE(J<=IVALUE)
            IF(J<IVALUE) THEN
              IF(A(KEY_DIMENSION,J)<A(KEY_DIMENSION,J+1)) J=J+1
            ENDIF
            IF(VALUE(KEY_DIMENSION)<A(KEY_DIMENSION,J)) THEN
              A(:,I)=A(:,J)
              I=J
              J=J+J
            ELSE
              J=IVALUE+1
            ENDIF
          ENDDO
          A(:,I)=VALUE
        ENDDO
      ENDIF
    ELSE
      localError="The specified key dimension of "//TRIM(NumberToVString(KEY_DIMENSION,"*",err,error))// &
        & " is invalid. The key dimension must be > 0 and <= "//TRIM(NumberToVString(SIZE(A,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    Call Exits("LIST_SORT_HEAP_C_INT2_ARRAY")
    RETURN
999 CALL Errors("LIST_SORT_HEAP_C_INT2_ARRAY",err,error)
    Call Exits("LIST_SORT_HEAP_C_INT2_ARRAY")
    RETURN 1
  END SUBROUTINE LIST_SORT_HEAP_C_INT2_ARRAY
  
  !
  !================================================================================================================================
  !

  !>Sorts a real single precision array of data dimension 1 list into assending order using the heap sort method.
  SUBROUTINE LIST_SORT_HEAP_SP1_ARRAY(A,err,error,*)
  
    !Argument variables
    REAL(SP), INTENT(INOUT) :: A(:) !<The list to sort
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: I,IVALUE,J,L
    REAL(SP) :: VALUE
    
    CALL Enters("LIST_SORT_HEAP_SP1_ARRAY",err,error,*999)

    IF(SIZE(A,1)>1) THEN      
      L=SIZE(A,1)/2+1
      IVALUE=SIZE(A,1)
      DO 
        IF(L>1) THEN
          L=L-1
          VALUE=A(L)
        ELSE
          VALUE=A(IVALUE)
          A(IVALUE)=A(1)
          IVALUE=IVALUE-1
          IF(IVALUE==1) THEN
            A(1)=VALUE
            EXIT
          ENDIF
        ENDIF
        I=L
        J=L+L
        DO WHILE(J<=IVALUE)
          IF(J<IVALUE) THEN
            IF(A(J)<A(J+1)) J=J+1
          ENDIF
          IF(VALUE<A(J)) THEN
            A(I)=A(J)
            I=J
            J=J+J
          ELSE
            J=IVALUE+1
          ENDIF
        ENDDO
        A(I)=VALUE
      ENDDO
    ENDIF

    Call Exits("LIST_SORT_HEAP_SP1_ARRAY")
    RETURN
999 CALL Errors("LIST_SORT_HEAP_SP1_ARRAY",err,error)
    Call Exits("LIST_SORT_HEAP_SP1_ARRAY")
    RETURN 1
  END SUBROUTINE LIST_SORT_HEAP_SP1_ARRAY
  
  !
  !================================================================================================================================
  !

  !>Sorts a real single precision array of data dimension > 1 list into assending order using the heap sort method.
  SUBROUTINE LIST_SORT_HEAP_SP2_ARRAY(A,KEY_DIMENSION,err,error,*)
  
    !Argument variables
    REAL(SP), INTENT(INOUT) :: A(:,:) !<The list to sort
    INTEGER(INTG), INTENT(IN) :: KEY_DIMENSION !<The key dimension of A to do the sort on
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: I,IVALUE,J,L
    REAL(SP) :: VALUE(SIZE(A,1))
    TYPE(VARYING_STRING) :: localError
    
    CALL Enters("LIST_SORT_HEAP_SP2_ARRAY",err,error,*999)

    IF(KEY_DIMENSION>0.AND.KEY_DIMENSION<=SIZE(A,1)) THEN
      IF(SIZE(A,2)>1) THEN      
        L=SIZE(A,2)/2+1
        IVALUE=SIZE(A,2)
        DO 
          IF(L>1) THEN
            L=L-1
            VALUE=A(:,L)
          ELSE
            VALUE=A(:,IVALUE)
            A(:,IVALUE)=A(:,1)
            IVALUE=IVALUE-1
            IF(IVALUE==1) THEN
              A(:,1)=VALUE
              EXIT
            ENDIF
          ENDIF
          I=L
          J=L+L
          DO WHILE(J<=IVALUE)
            IF(J<IVALUE) THEN
              IF(A(KEY_DIMENSION,J)<A(KEY_DIMENSION,J+1)) J=J+1
            ENDIF
            IF(VALUE(KEY_DIMENSION)<A(KEY_DIMENSION,J)) THEN
              A(:,I)=A(:,J)
              I=J
              J=J+J
            ELSE
              J=IVALUE+1
            ENDIF
          ENDDO
          A(:,I)=VALUE
        ENDDO
      ENDIF
    ELSE
      localError="The specified key dimension of "//TRIM(NumberToVString(KEY_DIMENSION,"*",err,error))// &
        & " is invalid. The key dimension must be > 0 and <= "//TRIM(NumberToVString(SIZE(A,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    Call Exits("LIST_SORT_HEAP_SP2_ARRAY")
    RETURN
999 CALL Errors("LIST_SORT_HEAP_SP2_ARRAY",err,error)
    Call Exits("LIST_SORT_HEAP_SP2_ARRAY")
    RETURN 1
  END SUBROUTINE LIST_SORT_HEAP_SP2_ARRAY
  
  !
  !================================================================================================================================
  !
  
  !>Sorts a real double precision array of data dimension 1 list into assending order using the heap sort method.
  SUBROUTINE LIST_SORT_HEAP_DP1_ARRAY(A,err,error,*)
  
    !Argument variables
    REAL(DP), INTENT(INOUT) :: A(:) !<The list to sort
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: I,IVALUE,J,L
    REAL(DP) :: VALUE
    
    CALL Enters("LIST_SORT_HEAP_DP1_ARRAY",err,error,*999)

    IF(SIZE(A,1)>1) THEN      
      L=SIZE(A,1)/2+1
      IVALUE=SIZE(A,1)
      DO 
        IF(L>1) THEN
          L=L-1
          VALUE=A(L)
        ELSE
          VALUE=A(IVALUE)
          A(IVALUE)=A(1)
          IVALUE=IVALUE-1
          IF(IVALUE==1) THEN
            A(1)=VALUE
            EXIT
          ENDIF
        ENDIF
        I=L
        J=L+L
        DO WHILE(J<=IVALUE)
          IF(J<IVALUE) THEN
            IF(A(J)<A(J+1)) J=J+1
          ENDIF
          IF(VALUE<A(J)) THEN
            A(I)=A(J)
            I=J
            J=J+J
          ELSE
            J=IVALUE+1
          ENDIF
        ENDDO
        A(I)=VALUE
      ENDDO
    ENDIF

    Call Exits("LIST_SORT_HEAP_DP1_ARRAY")
    RETURN
999 CALL Errors("LIST_SORT_HEAP_DP1_ARRAY",err,error)
    Call Exits("LIST_SORT_HEAP_DP1_ARRAY")
    RETURN 1
  END SUBROUTINE LIST_SORT_HEAP_DP1_ARRAY
  
  !
  !================================================================================================================================
  !
  
  !>Sorts a real double precision array of data dimension > 1 list into assending order using the heap sort method.
  SUBROUTINE LIST_SORT_HEAP_DP2_ARRAY(A,KEY_DIMENSION,err,error,*)
  
    !Argument variables
    REAL(DP), INTENT(INOUT) :: A(:,:) !<The list to sort
    INTEGER(INTG), INTENT(IN) :: KEY_DIMENSION !<The key dimension of A to do the sort on
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: I,IVALUE,J,L
    REAL(DP) :: VALUE(SIZE(A,1))
    TYPE(VARYING_STRING) :: localError
   
    CALL Enters("LIST_SORT_HEAP_DP2_ARRAY",err,error,*999)

    IF(KEY_DIMENSION>0.AND.KEY_DIMENSION<=SIZE(A,1)) THEN
      IF(SIZE(A,2)>1) THEN      
        L=SIZE(A,2)/2+1
        IVALUE=SIZE(A,2)
        DO 
          IF(L>1) THEN
            L=L-1
            VALUE=A(:,L)
          ELSE
            VALUE=A(:,IVALUE)
            A(:,IVALUE)=A(:,1)
            IVALUE=IVALUE-1
            IF(IVALUE==1) THEN
              A(:,1)=VALUE
              EXIT
            ENDIF
          ENDIF
          I=L
          J=L+L
          DO WHILE(J<=IVALUE)
            IF(J<IVALUE) THEN
              IF(A(KEY_DIMENSION,J)<A(KEY_DIMENSION,J+1)) J=J+1
            ENDIF
            IF(VALUE(KEY_DIMENSION)<A(KEY_DIMENSION,J)) THEN
              A(:,I)=A(:,J)
              I=J
              J=J+J
            ELSE
              J=IVALUE+1
            ENDIF
          ENDDO
          A(:,I)=VALUE
        ENDDO
      ENDIF
    ELSE
      localError="The specified key dimension of "//TRIM(NumberToVString(KEY_DIMENSION,"*",err,error))// &
        & " is invalid. The key dimension must be > 0 and <= "//TRIM(NumberToVString(SIZE(A,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    Call Exits("LIST_SORT_HEAP_DP2_ARRAY")
    RETURN
999 CALL Errors("LIST_SORT_HEAP_DP2_ARRAY",err,error)
    Call Exits("LIST_SORT_HEAP_DP2_ARRAY")
    RETURN 1
  END SUBROUTINE LIST_SORT_HEAP_DP2_ARRAY
  
  !
  !================================================================================================================================
  !

  !>Sorts an integer array of data dimension 1 list into either assending or descending order using the shell sort method.
  SUBROUTINE LIST_SORT_SHELL_INTG1_ARRAY(A,err,error,*)
  
    !Argument variables
    INTEGER(INTG), INTENT(INOUT) :: A(:) !<The list to sort
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: I,INC,J,VALUE
    
    CALL Enters("LIST_SORT_SHELL_INTG1_ARRAY",err,error,*999)

    INC=4
    DO WHILE(INC<=SIZE(A,1))
      INC=3*INC+1
    ENDDO
    DO WHILE(INC>1)
      INC=INC/3
      DO i=INC+1,SIZE(A,1)
        VALUE=A(i)
        J=I
        DO WHILE(A(J-INC)>VALUE)
          A(J)=A(J-INC)
          J=J-INC
          IF(J<=INC) EXIT
        ENDDO
        A(J)=VALUE
      ENDDO !i
    ENDDO

    Call Exits("LIST_SORT_SHELL_INTG1_ARRAY")
    RETURN
999 CALL Errors("LIST_SORT_SHELL_INTG1_ARRAY",err,error)
    Call Exits("LIST_SORT_SHELL_INTG1_ARRAY")
    RETURN 1
  END SUBROUTINE LIST_SORT_SHELL_INTG1_ARRAY
  
  !
  !================================================================================================================================
  !

  !>Sorts an integer array of data dimension > 1 list into either assending or descending order using the shell sort method.
  SUBROUTINE LIST_SORT_SHELL_INTG2_ARRAY(A,KEY_DIMENSION,err,error,*)
  
    !Argument variables
    INTEGER(INTG), INTENT(INOUT) :: A(:,:) !<The list to sort
    INTEGER(INTG), INTENT(IN) :: KEY_DIMENSION !<The key dimension of A to do the sort on
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: I,INC,J,VALUE(SIZE(A,1))
    TYPE(VARYING_STRING) :: localError
   
    CALL Enters("LIST_SORT_SHELL_INTG2_ARRAY",err,error,*999)

    IF(KEY_DIMENSION>0.AND.KEY_DIMENSION<=SIZE(A,1)) THEN
      INC=4
      DO WHILE(INC<=SIZE(A,2))
        INC=3*INC+1
      ENDDO
      DO WHILE(INC>1)
        INC=INC/3
        DO i=INC+1,SIZE(A,2)
          VALUE=A(:,i)
          J=I
          DO WHILE(A(KEY_DIMENSION,J-INC)>VALUE(KEY_DIMENSION))
            A(:,J)=A(:,J-INC)
            J=J-INC
            IF(J<=INC) EXIT
          ENDDO
          A(:,J)=VALUE
        ENDDO !i
      ENDDO
    ELSE
      localError="The specified key dimension of "//TRIM(NumberToVString(KEY_DIMENSION,"*",err,error))// &
        & " is invalid. The key dimension must be > 0 and <= "//TRIM(NumberToVString(SIZE(A,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    Call Exits("LIST_SORT_SHELL_INTG2_ARRAY")
    RETURN
999 CALL Errors("LIST_SORT_SHELL_INTG2_ARRAY",err,error)
    Call Exits("LIST_SORT_SHELL_INTG2_ARRAY")
    RETURN 1
  END SUBROUTINE LIST_SORT_SHELL_INTG2_ARRAY
  
  !
  !================================================================================================================================
  !

  !>Sorts an integer array of data dimension 1 list into either assending or descending order using the shell sort method.
  SUBROUTINE LIST_SORT_SHELL_C_INT1_ARRAY(A,err,error,*)
  
    !Argument variables
    INTEGER(C_INT), INTENT(INOUT) :: A(:) !<The list to sort
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: I,INC,J
    INTEGER(C_INT) :: VALUE
    
    CALL Enters("LIST_SORT_SHELL_C_INT1_ARRAY",err,error,*999)

    INC=4
    DO WHILE(INC<=SIZE(A,1))
      INC=3*INC+1
    ENDDO
    DO WHILE(INC>1)
      INC=INC/3
      DO i=INC+1,SIZE(A,1)
        VALUE=A(i)
        J=I
        DO WHILE(A(J-INC)>VALUE)
          A(J)=A(J-INC)
          J=J-INC
          IF(J<=INC) EXIT
        ENDDO
        A(J)=VALUE
      ENDDO !i
    ENDDO

    Call Exits("LIST_SORT_SHELL_C_INT1_ARRAY")
    RETURN
999 CALL Errors("LIST_SORT_SHELL_C_INT1_ARRAY",err,error)
    Call Exits("LIST_SORT_SHELL_C_INT1_ARRAY")
    RETURN 1
  END SUBROUTINE LIST_SORT_SHELL_C_INT1_ARRAY
  
  !
  !================================================================================================================================
  !

  !>Sorts an integer array of data dimension > 1 list into either assending or descending order using the shell sort method.
  SUBROUTINE LIST_SORT_SHELL_C_INT2_ARRAY(A,KEY_DIMENSION,err,error,*)
  
    !Argument variables
    INTEGER(C_INT), INTENT(INOUT) :: A(:,:) !<The list to sort
    INTEGER(INTG), INTENT(IN) :: KEY_DIMENSION !<The key dimension of A to do the sort on
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: I,INC,J
    INTEGER(C_INT) :: VALUE(SIZE(A,1))
    TYPE(VARYING_STRING) :: localError
   
    CALL Enters("LIST_SORT_SHELL_C_INT2_ARRAY",err,error,*999)

    IF(KEY_DIMENSION>0.AND.KEY_DIMENSION<=SIZE(A,1)) THEN
      INC=4
      DO WHILE(INC<=SIZE(A,2))
        INC=3*INC+1
      ENDDO
      DO WHILE(INC>1)
        INC=INC/3
        DO i=INC+1,SIZE(A,2)
          VALUE=A(:,i)
          J=I
          DO WHILE(A(KEY_DIMENSION,J-INC)>VALUE(KEY_DIMENSION))
            A(:,J)=A(:,J-INC)
            J=J-INC
            IF(J<=INC) EXIT
          ENDDO
          A(:,J)=VALUE
        ENDDO !i
      ENDDO
    ELSE
      localError="The specified key dimension of "//TRIM(NumberToVString(KEY_DIMENSION,"*",err,error))// &
        & " is invalid. The key dimension must be > 0 and <= "//TRIM(NumberToVString(SIZE(A,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    Call Exits("LIST_SORT_SHELL_C_INT2_ARRAY")
    RETURN
999 CALL Errors("LIST_SORT_SHELL_C_INT2_ARRAY",err,error)
    Call Exits("LIST_SORT_SHELL_C_INT2_ARRAY")
    RETURN 1
  END SUBROUTINE LIST_SORT_SHELL_C_INT2_ARRAY
  
  !
  !================================================================================================================================
  !

  !>Sorts a real single precision array of data dimension 1 list into either assending or descending order using the shell
  !>sort method.
  SUBROUTINE LIST_SORT_SHELL_SP1_ARRAY(A,err,error,*)
  
    !Argument variables
    REAL(SP), INTENT(INOUT) :: A(:) !<The list to sort
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: I,INC,J
    REAL(SP) :: VALUE
    
    CALL Enters("LIST_SORT_SHELL_SP1_ARRAY",err,error,*999)

    INC=4
    DO WHILE(INC<=SIZE(A,1))
      INC=3*INC+1
    ENDDO
    DO WHILE(INC>1)
      INC=INC/3
      DO i=INC+1,SIZE(A,1)
        VALUE=A(i)
        J=I
        DO WHILE(A(J-INC)>VALUE)
          A(J)=A(J-INC)
          J=J-INC
          IF(J<=INC) EXIT
        ENDDO
        A(J)=VALUE
      ENDDO !i
    ENDDO

    Call Exits("LIST_SORT_SHELL_SP1_ARRAY")
    RETURN
999 CALL Errors("LIST_SORT_SHELL_SP1_ARRAY",err,error)
    Call Exits("LIST_SORT_SHELL_SP1_ARRAY")
    RETURN 1
  END SUBROUTINE LIST_SORT_SHELL_SP1_ARRAY
  
  !
  !================================================================================================================================
  !

  !>Sorts a real single precision array of data dimension > 1 list into either assending or descending order using the shell
  !>sort method.
  SUBROUTINE LIST_SORT_SHELL_SP2_ARRAY(A,KEY_DIMENSION,err,error,*)
  
    !Argument variables
    REAL(SP), INTENT(INOUT) :: A(:,:) !<The list to sort
    INTEGER(INTG), INTENT(IN) :: KEY_DIMENSION !<The key dimension of A to do the sort on
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: I,INC,J
    REAL(SP) :: VALUE(SIZE(A,1))
    TYPE(VARYING_STRING) :: localError
    
    CALL Enters("LIST_SORT_SHELL_SP2_ARRAY",err,error,*999)

    IF(KEY_DIMENSION>0.AND.KEY_DIMENSION<=SIZE(A,1)) THEN
      INC=4
      DO WHILE(INC<=SIZE(A,2))
        INC=3*INC+1
      ENDDO
      DO WHILE(INC>1)
        INC=INC/3
        DO i=INC+1,SIZE(A,2)
          VALUE=A(:,i)
          J=I
          DO WHILE(A(KEY_DIMENSION,J-INC)>VALUE(KEY_DIMENSION))
            A(:,J)=A(:,J-INC)
            J=J-INC
            IF(J<=INC) EXIT
          ENDDO
          A(:,J)=VALUE
        ENDDO !i
      ENDDO
    ELSE
      localError="The specified key dimension of "//TRIM(NumberToVString(KEY_DIMENSION,"*",err,error))// &
        & " is invalid. The key dimension must be > 0 and <= "//TRIM(NumberToVString(SIZE(A,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    Call Exits("LIST_SORT_SHELL_SP2_ARRAY")
    RETURN
999 CALL Errors("LIST_SORT_SHELL_SP2_ARRAY",err,error)
    Call Exits("LIST_SORT_SHELL_SP2_ARRAY")
    RETURN 1
  END SUBROUTINE LIST_SORT_SHELL_SP2_ARRAY
  
  !
  !================================================================================================================================
  !

  !>Sorts a real double precision array of data dimension 1 list into either assending or descending order using the shell
  !>sort method.
  SUBROUTINE LIST_SORT_SHELL_DP1_ARRAY(A,err,error,*)
  
    !Argument variables
    REAL(DP), INTENT(INOUT) :: A(:)
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: I,INC,J
    REAL(DP) :: VALUE
    
    CALL Enters("LIST_SORT_SHELL_DP1_ARRAY",err,error,*999)

    INC=4
    DO WHILE(INC<=SIZE(A,1))
      INC=3*INC+1
    ENDDO
    DO WHILE(INC>1)
      INC=INC/3
      DO i=INC+1,SIZE(A,1)
        VALUE=A(i)
        J=I
        DO WHILE(A(J-INC)>VALUE)
          A(J)=A(J-INC)
          J=J-INC
          IF(J<=INC) EXIT
        ENDDO
        A(J)=VALUE
      ENDDO !i
    ENDDO

    Call Exits("LIST_SORT_SHELL_DP1_ARRAY")
    RETURN
999 CALL Errors("LIST_SORT_SHELL_DP1_ARRAY",err,error)
    Call Exits("LIST_SORT_SHELL_DP1_ARRAY")
    RETURN 1
  END SUBROUTINE LIST_SORT_SHELL_DP1_ARRAY
  
  !
  !================================================================================================================================
  !

  !>Sorts a real double precision array of data dimension 2 list into either assending or descending order using the shell
  !>sort method.
  SUBROUTINE LIST_SORT_SHELL_DP2_ARRAY(A,KEY_DIMENSION,err,error,*)
  
    !Argument variables
    REAL(DP), INTENT(INOUT) :: A(:,:) !<The list to sort
    INTEGER(INTG), INTENT(IN) :: KEY_DIMENSION !<The key dimension of A to do the sort on
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: I,INC,J
    REAL(DP) :: VALUE(SIZE(A,1))
    TYPE(VARYING_STRING) :: localError
    
    CALL Enters("LIST_SORT_SHELL_DP2_ARRAY",err,error,*999)

    IF(KEY_DIMENSION>0.AND.KEY_DIMENSION<=SIZE(A,1)) THEN
      INC=4
      DO WHILE(INC<=SIZE(A,2))
        INC=3*INC+1
      ENDDO
      DO WHILE(INC>1)
        INC=INC/3
        DO i=INC+1,SIZE(A,2)
          VALUE=A(:,i)
          J=I
          DO WHILE(A(KEY_DIMENSION,J-INC)>VALUE(KEY_DIMENSION))
            A(:,J)=A(:,J-INC)
            J=J-INC
            IF(J<=INC) EXIT
          ENDDO
          A(:,J)=VALUE
        ENDDO !i
      ENDDO
    ELSE
      localError="The specified key dimension of "//TRIM(NumberToVString(KEY_DIMENSION,"*",err,error))// &
        & " is invalid. The key dimension must be > 0 and <= "//TRIM(NumberToVString(SIZE(A,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    Call Exits("LIST_SORT_SHELL_DP1_ARRAY")
    RETURN
999 CALL Errors("LIST_SORT_SHELL_DP1_ARRAY",err,error)
    Call Exits("LIST_SORT_SHELL_DP1_ARRAY")
    RETURN 1
  END SUBROUTINE LIST_SORT_SHELL_DP2_ARRAY
  
  !
  !================================================================================================================================
  !

  !>Finds the intersection of two sets (arrays), leaving the original arrays intact
  SUBROUTINE LIST_INTERSECTION_INTG_ARRAY(A,B,C,err,error,*)
    
    ! Argument variables
    INTEGER(INTG), INTENT(IN), TARGET :: A(:)   !<One of the two arrays to find the intersection for
    INTEGER(INTG), INTENT(IN), TARGET :: B(:)   !<Other array to find the intersection for
    INTEGER(INTG), ALLOCATABLE, INTENT(OUT) :: C(:) !<On exit, contains the list of common elements of the arrays
    INTEGER(INTG), INTENT(OUT) :: ERR          !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    ! Local variables
    INTEGER(INTG) :: SIZE_SHORTER,SIZE_LONGER
    INTEGER(INTG) :: I,J,START,NUMBER_OF_MATCHES
    INTEGER(INTG), POINTER :: LONGER(:),SHORTER(:)
    INTEGER(INTG), ALLOCATABLE :: MATCHES(:)
    INTEGER(INTG), ALLOCATABLE :: LONG_ARRAY(:),SHORT_ARRAY(:)   !<copies, if needed
    
    CALL Enters("LIST_INTERSECTION_INTG_ARRAY",err,error,*999)

    ! if the lists are small, it's probably easier to directly compare: O(n^2)
    ! but if they're big, sort first then compare: O(n log n)*2 + 2*O(n)

    IF(ALLOCATED(C)) THEN
      ! theoretically this cannot happen?
      CALL FlagError("Output array is already allocated.",err,error,*999)
    ELSE
      ! start finding the intersection
      NULLIFY(LONGER)
      NULLIFY(SHORTER)
      ! it's quicker to compare shorter array elements to longer ones
      IF(SIZE(A)>SIZE(B)) THEN
        LONGER=>A
        SHORTER=>B
      ELSE
        LONGER=>B
        SHORTER=>A
      ENDIF
      SIZE_SHORTER=SIZE(SHORTER)
      SIZE_LONGER=SIZE(LONGER)
      ALLOCATE(MATCHES(SIZE_SHORTER))
      NUMBER_OF_MATCHES=0

      ! long or short lists?
      IF(SIZE_LONGER*SIZE_SHORTER<=1E4) THEN  ! a rather arbitrary cutoff...
        ! 'short' lists - begin comparing straight away
        DO I=1,SIZE_SHORTER
          DO J=1,SIZE_LONGER
            IF(SHORTER(I)==LONGER(J)) THEN
              NUMBER_OF_MATCHES=NUMBER_OF_MATCHES+1
              MATCHES(NUMBER_OF_MATCHES)=SHORTER(I)
            ENDIF
          ENDDO
        ENDDO
      ELSE
        ! 'long' lists - make copies of the arrays
        ALLOCATE(LONG_ARRAY(SIZE_LONGER),SHORT_ARRAY(SIZE_SHORTER))
        LONG_ARRAY=LONGER
        SHORT_ARRAY=SHORTER
        ! sort both arrays
        CALL LIST_SORT(LONG_ARRAY,err,error,*999)
        CALL LIST_SORT(SHORT_ARRAY,err,error,*999)
        ! compare now
        START=1
        DO I=1,SIZE_SHORTER
          DO J=START,SIZE_LONGER
            IF(LONG_ARRAY(J)==SHORT_ARRAY(I)) THEN
              NUMBER_OF_MATCHES=NUMBER_OF_MATCHES+1
              MATCHES(NUMBER_OF_MATCHES)=SHORT_ARRAY(I)
              START=MIN(J+1,SIZE_LONGER)
              EXIT
            ELSEIF(LONG_ARRAY(J)>SHORT_ARRAY(I)) THEN
              ! can start here next time
              START=MAX(J-1,1)
              EXIT
            ENDIF
          ENDDO
        ENDDO
        DEALLOCATE(LONG_ARRAY,SHORT_ARRAY)
      ENDIF ! long or short lists
      ! cut the array down to size
      ALLOCATE(C(NUMBER_OF_MATCHES))
      C=MATCHES(1:NUMBER_OF_MATCHES)
      DEALLOCATE(MATCHES)
    ENDIF

    Call Exits("LIST_INTERSECTION_INTG_ARRAY")
    RETURN
999 CALL Errors("LIST_INTERSECTION_INTG_ARRAY",err,error)
    Call Exits("LIST_INTERSECTION_INTG_ARRAY")
    RETURN 1

  END SUBROUTINE LIST_INTERSECTION_INTG_ARRAY

  !
  !================================================================================================================================
  !

  !>Finds the intersection of two sets (arrays), leaving the original arrays intact
  SUBROUTINE LIST_INTERSECTION_C_INT_ARRAY(A,B,C,err,error,*)
    
    ! Argument variables
    INTEGER(C_INT), INTENT(IN), TARGET :: A(:)   !<One of the two arrays to find the intersection for
    INTEGER(C_INT), INTENT(IN), TARGET :: B(:)   !<Other array to find the intersection for
    INTEGER(C_INT), ALLOCATABLE, INTENT(OUT) :: C(:) !<On exit, contains the list of common elements of the arrays
    INTEGER(INTG), INTENT(OUT) :: ERR          !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    ! Local variables
    INTEGER(INTG) :: SIZE_SHORTER,SIZE_LONGER
    INTEGER(INTG) :: I,J,START,NUMBER_OF_MATCHES
    INTEGER(C_INT), POINTER :: LONGER(:),SHORTER(:)
    INTEGER(C_INT), ALLOCATABLE :: MATCHES(:)
    INTEGER(C_INT), ALLOCATABLE :: LONG_ARRAY(:),SHORT_ARRAY(:)   !<copies, if needed
    
    CALL Enters("LIST_INTERSECTION_C_INT_ARRAY",err,error,*999)

    ! if the lists are small, it's probably easier to directly compare: O(n^2)
    ! but if they're big, sort first then compare: O(n log n)*2 + 2*O(n)

    IF(ALLOCATED(C)) THEN
      ! theoretically this cannot happen?
      CALL FlagError("Output array is already allocated.",err,error,*999)
    ELSE
      ! start finding the intersection
      NULLIFY(LONGER)
      NULLIFY(SHORTER)
      ! it's quicker to compare shorter array elements to longer ones
      IF(SIZE(A)>SIZE(B)) THEN
        LONGER=>A
        SHORTER=>B
      ELSE
        LONGER=>B
        SHORTER=>A
      ENDIF
      SIZE_SHORTER=SIZE(SHORTER)
      SIZE_LONGER=SIZE(LONGER)
      ALLOCATE(MATCHES(SIZE_SHORTER))
      NUMBER_OF_MATCHES=0

      ! long or short lists?
      IF(SIZE_LONGER*SIZE_SHORTER<=1E4) THEN  ! a rather arbitrary cutoff...
        ! 'short' lists - begin comparing straight away
        DO I=1,SIZE_SHORTER
          DO J=1,SIZE_LONGER
            IF(SHORTER(I)==LONGER(J)) THEN
              NUMBER_OF_MATCHES=NUMBER_OF_MATCHES+1
              MATCHES(NUMBER_OF_MATCHES)=SHORTER(I)
            ENDIF
          ENDDO
        ENDDO
      ELSE
        ! 'long' lists - make copies of the arrays
        ALLOCATE(LONG_ARRAY(SIZE_LONGER),SHORT_ARRAY(SIZE_SHORTER))
        LONG_ARRAY=LONGER
        SHORT_ARRAY=SHORTER
        ! sort both arrays
        CALL LIST_SORT(LONG_ARRAY,err,error,*999)
        CALL LIST_SORT(SHORT_ARRAY,err,error,*999)
        ! compare now
        START=1
        DO I=1,SIZE_SHORTER
          DO J=START,SIZE_LONGER
            IF(LONG_ARRAY(J)==SHORT_ARRAY(I)) THEN
              NUMBER_OF_MATCHES=NUMBER_OF_MATCHES+1
              MATCHES(NUMBER_OF_MATCHES)=SHORT_ARRAY(I)
              START=MIN(J+1,SIZE_LONGER)
              EXIT
            ELSEIF(LONG_ARRAY(J)>SHORT_ARRAY(I)) THEN
              ! can start here next time
              START=MAX(J-1,1)
              EXIT
            ENDIF
          ENDDO
        ENDDO
        DEALLOCATE(LONG_ARRAY,SHORT_ARRAY)
      ENDIF ! long or short lists
      ! cut the array down to size
      ALLOCATE(C(NUMBER_OF_MATCHES))
      C=MATCHES(1:NUMBER_OF_MATCHES)
      DEALLOCATE(MATCHES)
    ENDIF

    Call Exits("LIST_INTERSECTION_C_INT_ARRAY")
    RETURN
999 CALL Errors("LIST_INTERSECTION_C_INT_ARRAY",err,error)
    Call Exits("LIST_INTERSECTION_C_INT_ARRAY")
    RETURN 1

  END SUBROUTINE LIST_INTERSECTION_C_INT_ARRAY

  !
  !================================================================================================================================
  !

  !>Finds out whether array A is a subset of array B
  SUBROUTINE LISTS_SUBSET_OF_INTG_ARRAY(A,B,SUBSET,err,error,*)
    ! Argument variables
    INTEGER(INTG), INTENT(IN) :: A(:)   !<Supposed subset (to test for)
    INTEGER(INTG), INTENT(IN) :: B(:)   !<Supposed superset
    LOGICAL, INTENT(OUT) :: SUBSET              !<On exit, TRUE if A is a subset of B
    INTEGER(INTG), INTENT(OUT) :: ERR           !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR  !<The error string
    ! Logical variables
    INTEGER(INTG) :: SIZE_A,SIZE_B,I,J,START,SIZE_REDUCE
    INTEGER(INTG), ALLOCATABLE :: A_SORTED(:),B_SORTED(:)

    CALL Enters("LISTS_SUBSET_OF_INTG_ARRAY",err,error,*999)

    SIZE_A=SIZE(A)
    SIZE_B=SIZE(B)
    SUBSET=.FALSE.
    
    ! some easy tests
    IF(SIZE_A>SIZE_B) THEN
      Call Exits("LISTS_SUBSET_OF_INTG_ARRAY")
      RETURN
    ENDIF

    SIZE_REDUCE=0
    DO I=1,SIZE_A
      IF(A(I)==0) SIZE_REDUCE=SIZE_REDUCE+1
    ENDDO
    SIZE_A=SIZE_A-SIZE_REDUCE
    SIZE_REDUCE=0
    DO I=1,SIZE_B
      IF(B(I)==0) SIZE_REDUCE=SIZE_REDUCE+1
    ENDDO
    SIZE_B=SIZE_B-SIZE_REDUCE

    ! short of long arrays?
    IF(SIZE_A*SIZE_B<=1E4) THEN
      ! 'short' arrays - just compare without sorting
      DO I=1,SIZE_A
        DO J=1,SIZE_B
          IF(A(I)==B(J)) THEN
            EXIT
          ELSEIF(J==SIZE_B) THEN
            Call Exits("LISTS_SUBSET_OF_INTG_ARRAY")
            RETURN
          ENDIF
        ENDDO
        IF(I==SIZE_A) SUBSET=.TRUE.
      ENDDO
    ELSE
      ! 'long' arrays - sort first
      ALLOCATE(A_SORTED(SIZE_A),B_SORTED(SIZE_B))
      A_SORTED=A
      B_SORTED=B
      CALL LIST_SORT(A_SORTED,err,error,*999)
      CALL LIST_SORT(B_SORTED,err,error,*999)
      START=1
      DO I=1,SIZE_A
        DO J=1,SIZE_B
          IF(A(I)==B(J)) THEN
            START=MIN(J+1,SIZE_B)
            EXIT
          ELSEIF(A(I)<B(J)) THEN
            DEALLOCATE(A_SORTED,B_SORTED)
            Call Exits("LISTS_SUBSET_OF_INTG_ARRAY")
            RETURN
          ENDIF
        ENDDO
        IF(I==SIZE_A) SUBSET=.TRUE.
      ENDDO
      DEALLOCATE(A_SORTED,B_SORTED)
    ENDIF

    Call Exits("LISTS_SUBSET_OF_INTG_ARRAY")
    RETURN
999 CALL Errors("LISTS_SUBSET_OF_INTG_ARRAY",err,error)
    Call Exits("LISTS_SUBSET_OF_INTG_ARRAY")
    RETURN 1

  END SUBROUTINE LISTS_SUBSET_OF_INTG_ARRAY

  !
  !================================================================================================================================
  !

  !>Finds out whether array A is a subset of array B
  SUBROUTINE LISTS_SUBSET_OF_C_INT_ARRAY(A,B,SUBSET,err,error,*)
    ! Argument variables
    INTEGER(C_INT), INTENT(IN) :: A(:)   !<Supposed subset (to test for)
    INTEGER(C_INT), INTENT(IN) :: B(:)   !<Supposed superset
    LOGICAL, INTENT(OUT) :: SUBSET              !<On exit, TRUE if A is a subset of B
    INTEGER(INTG), INTENT(OUT) :: ERR           !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR  !<The error string
    ! Logical variables
    INTEGER(INTG) :: SIZE_A,SIZE_B,I,J,START,SIZE_REDUCE
    INTEGER(C_INT), ALLOCATABLE :: A_SORTED(:),B_SORTED(:)

    CALL Enters("LISTS_SUBSET_OF_C_INT_ARRAY",err,error,*999)

    SIZE_A=SIZE(A)
    SIZE_B=SIZE(B)
    SUBSET=.FALSE.
    
    ! some easy tests
    IF(SIZE_A>SIZE_B) THEN
      Call Exits("LISTS_SUBSET_OF_C_INT_ARRAY")
      RETURN
    ENDIF

    SIZE_REDUCE=0
    DO I=1,SIZE_A
      IF(A(I)==0) SIZE_REDUCE=SIZE_REDUCE+1
    ENDDO
    SIZE_A=SIZE_A-SIZE_REDUCE
    SIZE_REDUCE=0
    DO I=1,SIZE_B
      IF(B(I)==0) SIZE_REDUCE=SIZE_REDUCE+1
    ENDDO
    SIZE_B=SIZE_B-SIZE_REDUCE

    ! short of long arrays?
    IF(SIZE_A*SIZE_B<=1E4) THEN
      ! 'short' arrays - just compare without sorting
      DO I=1,SIZE_A
        DO J=1,SIZE_B
          IF(A(I)==B(J)) THEN
            EXIT
          ELSEIF(J==SIZE_B) THEN
            Call Exits("LISTS_SUBSET_OF_C_INT_ARRAY")
            RETURN
          ENDIF
        ENDDO
        IF(I==SIZE_A) SUBSET=.TRUE.
      ENDDO
    ELSE
      ! 'long' arrays - sort first
      ALLOCATE(A_SORTED(SIZE_A),B_SORTED(SIZE_B))
      A_SORTED=A
      B_SORTED=B
      CALL LIST_SORT(A_SORTED,err,error,*999)
      CALL LIST_SORT(B_SORTED,err,error,*999)
      START=1
      DO I=1,SIZE_A
        DO J=1,SIZE_B
          IF(A(I)==B(J)) THEN
            START=MIN(J+1,SIZE_B)
            EXIT
          ELSEIF(A(I)<B(J)) THEN
            DEALLOCATE(A_SORTED,B_SORTED)
            Call Exits("LISTS_SUBSET_OF_C_INT_ARRAY")
            RETURN
          ENDIF
        ENDDO
        IF(I==SIZE_A) SUBSET=.TRUE.
      ENDDO
      DEALLOCATE(A_SORTED,B_SORTED)
    ENDIF

    Call Exits("LISTS_SUBSET_OF_C_INT_ARRAY")
    RETURN
999 CALL Errors("LISTS_SUBSET_OF_C_INT_ARRAY",err,error)
    Call Exits("LISTS_SUBSET_OF_C_INT_ARRAY")
    RETURN 1

  END SUBROUTINE LISTS_SUBSET_OF_C_INT_ARRAY

  !
  !================================================================================================================================
  !

END MODULE Lists
