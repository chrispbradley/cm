!> \file
!> \author Chris Bradley
!> \brief This module handles all solver mapping routines.
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

!> This module handles all solver mapping routines.
MODULE SolverMappingRoutines

  USE BASE_ROUTINES
  USE BOUNDARY_CONDITIONS_ROUTINES
  USE COMP_ENVIRONMENT
  USE DISTRIBUTED_MATRIX_VECTOR
  USE DOMAIN_MAPPINGS
  USE EQUATIONS_SET_CONSTANTS
  USE FIELD_ROUTINES
  USE INPUT_OUTPUT
  USE INTERFACE_CONDITIONS_CONSTANTS
  USE ISO_VARYING_STRING
  USE KINDS
  USE LISTS
  USE MATRIX_VECTOR
  USE PROBLEM_CONSTANTS
  USE STRINGS
  USE TYPES

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !> \addtogroup SolverMapping_EquationsMatrixTypes SolverMapping::EquationsMatrixTypes
  !> \brief Equations matrix types
  !> \see SolverMapping
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_MAPPING_EQUATIONS_DYNAMIC_MATRIX=1 !<The equations matrix in the solver mapping is a dynamic equations matrix \see SolverMapping_EquationsMatrixTypes,SolverMapping
  INTEGER(INTG), PARAMETER :: SOLVER_MAPPING_EQUATIONS_LINEAR_MATRIX=2 !<The equations matrix in the solver mapping is a linear equations matrix \see SolverMapping_EquationsMatrixTypes,SolverMapping
  INTEGER(INTG), PARAMETER :: SOLVER_MAPPING_EQUATIONS_NONLINEAR_MATRIX=3 !<The equations matrix in the solver mapping is a nonlinear equations (Jacobian) matrix \see SolverMapping_EquationsMatrixTypes,SolverMapping
  !>@}
 
  !> \addtogroup SolverMapping_EquationsTypes SolverMapping::EquationsTypes
  !> \brief Equations Matrix types
  !> \see SolverMapping
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_MAPPING_EQUATIONS_EQUATIONS_SET=1 !<The equations in the solver mapping is from an equations set \see SolverMapping_EquationsTypes,SolverMapping
  INTEGER(INTG), PARAMETER :: SOLVER_MAPPING_EQUATIONS_INTERFACE_CONDITION=2 !<The equations in the solver mapping is from an interface condition \see SolverMapping_EquationsTypes,SolverMapping
  INTEGER(INTG), PARAMETER :: SOLVER_MAPPING_EQUATIONS_INTERFACE_TRANSPOSE=3 !<The equations in the solver mapping is from a transposed interface condition \see SolverMapping_EquationsTypes,SolverMapping
  !>@}
 
  !Module types

  !Module variables

  !Interfaces

  INTERFACE SolverMapping_EquationsVariablesToSolverMatrixSet
    MODULE PROCEDURE SolverMapping_EquationsVariablesToSolverMatrixSet0
    MODULE PROCEDURE SolverMapping_EquationsVariablesToSolverMatrixSet1
  END INTERFACE SolverMapping_EquationsVariablesToSolverMatrixSet
  
  PUBLIC SOLVER_MAPPING_EQUATIONS_DYNAMIC_MATRIX,SOLVER_MAPPING_EQUATIONS_LINEAR_MATRIX

  PUBLIC SOLVER_MAPPING_EQUATIONS_EQUATIONS_SET,SOLVER_MAPPING_EQUATIONS_INTERFACE_CONDITION

  PUBLIC SolverMapping_CreateFinish,SolverMapping_CreateStart

  PUBLIC SolverMapping_Destroy
  
  PUBLIC SolverMapping_EquationsSetAdd

  PUBLIC SolverMapping_InterfaceConditionAdd

  PUBLIC SolverMapping_EquationsVariablesToSolverMatrixSet

  PUBLIC SolverMapping_SolverMatricesNumberSet
  
CONTAINS

  !
  !=================================================================================================================================
  !
  
  !>Calculates the solver mappings
  SUBROUTINE SolverMapping_Calculate(solverMapping,err,error,*)

    !Argument variables
    TYPE(SolverMappingType), POINTER :: solverMapping !<A pointer to the solver mapping to calcualte
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: colEquationsColIdx,columnIdx,columnListItem(5),columnRank,dofIdx,dofType,dummyErr,dynamicVariableType, &
      & eqnLocalDof,equationType,equationsColumn,equationsIdx,equationsIdx2,equationsMatrix,equationsMatrixIdx,equationsRow, &
      & equationsRowNumber,equationsSetIdx,equationsVariableListItem(3),globalColumn,globalDof,globalDofCouplingNumber, &
      & globalDofIdx,globalDofsOffset,globalRow,globalRowIdx,interfaceColumn,interfaceColNumber,interfaceConditionIdx, &
      & interfaceConditionIdx2,interfaceEquationsListItem(2),interfaceIdx,interfaceMatrixIdx,interfaceRow,interfaceRowNumber, &
      & jacobianColumn,linearVariableType,lhsVariableType,localColumn,localDof,localDofsOffset,localRow,matricesType,matrixNumber, &
      & matrixType,matrixTypeIdx,matrixVariableIdx,myRank,numberOfColumns,numberOfDynamicEquationsMatrices, &
      & numberOfEquationsColumns,numberOfEquationsSets,numberOfEquationsVariables,numberofInterfaces,numberOfInterfaceColumns, &
      & numberOfInterfaceRows,numberOfInterfaceVariables,numberOfGlobalSolverDofs,numberOfGlobalSolverRows, &
      & numberOfLinearEquationsMatrices,numberOfLocalSolverDofs,numberOfLocalSolverRows,numberOfRankCols, &
      & numberOfRankRows,numberOfVariables,numberColEquationsCols,numberRowEquationsRows,rank,rankIdx,residualVariableType, &
      & rowEquationsRowIdx,rowIdx,rowListItem(4),rowRank,solverGlobalDof,solverMatrixIdx,solverVariableIdx,solverVariableIdxTemp, &
      & tempOffset,totalNumberOfLocalSolverDofs,variableIdx,variableListItem(3),variablePositionIdx,variableType
    INTEGER(INTG), ALLOCATABLE :: equationsSetVariables(:,:),equationsVariables(:,:),interfaceEquationsList(:,:), &
      & interfaceVariables(:,:),rankGlobalRowsList(:,:),rankGlobalColsList(:,:),solverLocalDof(:)
    INTEGER(INTG), ALLOCATABLE :: numberOfVariableGlobalSolverDofs(:),numberOfVariableLocalSolverDofs(:), &
      & totalNumberOfVariableLocalSolverDofs(:),subMatrixInformation(:,:,:),subMatrixList(:,:,:),variableTypes(:)
    LOGICAL :: constrainedDof,found,includeColumn,includeRow
    LOGICAL, ALLOCATABLE :: variableProcessed(:),variableRankProcessed(:,:)
    REAL(DP) :: couplingCoefficient
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: boundaryConditions
    TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER :: boundaryConditionsVariable
    TYPE(BoundaryConditionsCoupledDofsType), POINTER :: colEquationCols,rowEquationRows
    TYPE(BoundaryConditionsCoupledDofsType), TARGET :: dummyDofCoupling
    TYPE(BoundaryConditionsDofConstraintsType), POINTER :: dofConstraints
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: colDomainMapping,colDofsMapping,rowDomainMapping,rowDofsMapping
    TYPE(EQUATIONS_TYPE), POINTER :: equations
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: equationsMapping
    TYPE(EQUATIONS_MAPPING_DYNAMIC_TYPE), POINTER :: dynamicMapping
    TYPE(EQUATIONS_MAPPING_LINEAR_TYPE), POINTER :: linearMapping
    TYPE(EQUATIONS_MAPPING_NONLINEAR_TYPE), POINTER :: nonlinearMapping
    TYPE(EquationsMapping_LhsType), POINTER :: lhsMapping
    TYPE(EQUATIONS_MAPPING_RHS_TYPE), POINTER :: rhsMapping
    TYPE(EQUATIONS_MAPPING_SOURCE_TYPE), POINTER :: sourceMapping
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet
    TYPE(EquationsToSolverMapsType), POINTER :: equationsToSolverMap
    TYPE(FIELD_TYPE), POINTER :: dependentField,LagrangeField
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: dependentVariable,dynamicVariable,lagrangeVariable,lhsVariable,linearVariable, &
      & residualVariable,variable
    TYPE(INTEGER_INTG_PTR_TYPE), POINTER :: dofMap(:)
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: interfaceCondition
    TYPE(INTERFACE_DEPENDENT_TYPE), POINTER :: interfaceDependent
    TYPE(INTERFACE_EQUATIONS_TYPE), POINTER :: interfaceEquations
    TYPE(INTERFACE_MAPPING_TYPE), POINTER :: interfaceMapping
    TYPE(InterfaceToSolverMapsType), POINTER :: interfaceToSolverMap
    TYPE(JacobianToSolverMapType), POINTER :: jacobianToSolverMap
    TYPE(LIST_TYPE), POINTER :: equationsSetVariableList
    TYPE(LIST_PTR_TYPE), ALLOCATABLE :: interfaceEquationsLists(:),rankGlobalRowsLists(:,:),rankGlobalColsLists(:,:,:,:), &
      & variablesList(:)
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: solverEquations
    TYPE(SolverMappingDofCouplingsType) :: columnCouplings,rowCouplings
    TYPE(VARYING_STRING) :: dummyError,localError

    CALL Enters("SolverMapping_Calculate",err,error,*999)

    IF(ASSOCIATED(solverMapping)) THEN
      IF(ASSOCIATED(solverMapping%createValuesCache)) THEN
        solverEquations=>solverMapping%solverEquations
        IF(ASSOCIATED(solverEquations)) THEN
          boundaryConditions=>solverEquations%BOUNDARY_CONDITIONS
          IF(.NOT.ASSOCIATED(boundaryConditions)) &
            & CALL FlagError("The solver equations boundary conditions are not associated.",err,error,*999)
          
          !
          !--- Equations set <-> interface conditions  ---
          !
          ! 1. Calculate the lists of variables that are mapped to the rows and columns. From this calculate the sub-matrix
          !    map i.e., the locations of all the sub-matrices in the solver matrix. This will also determine which interface
          !    conditions influence an equations set and vice versa.
          !
          
          !Allocate equations set to solver map
          ALLOCATE(solverMapping%equationsSetToSolverMap(solverMapping%numberOfEquationsSets),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate solver mapping equations set to solver map.",err,error,*999)      
          !Allocate interface condition to solver map
          ALLOCATE(solverMapping%interfaceConditionToSolverMap(solverMapping%numberOfInterfaceConditions),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate solver mapping interface condition to solver map.",err,error,*999)
          !Allocate the column variable lists
          ALLOCATE(solverMapping%columnVariablesList(solverMapping%numberOfSolverMatrices),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate solver mapping column variables list.",err,error,*999)
          CALL SolverMapping_VariablesInitialise(solverMapping%rowVariablesList,err,error,*999)
          DO solverMatrixIdx=1,solverMapping%numberOfSolverMatrices
            CALL SolverMapping_VariablesInitialise(solverMapping%columnVariablesList(solverMatrixIdx),err,error,*999)
          ENDDO !solverMatrixIdx
          
          DO equationsSetIdx=1,solverMapping%numberOfEquationsSets
            equationsSet=>solverMapping%equationsSets(equationsSetIdx)%ptr
            IF(ASSOCIATED(equationsSet)) THEN
              equations=>equationsSet%equations
              IF(ASSOCIATED(equations)) THEN
                CALL SolverMapping_EquationsSetToSolverMapInitialise(solverMapping%equationsSetToSolverMap(equationsSetIdx), &
                  & err,error,*999)
                solverMapping%equationsSetToSolverMap(equationsSetIdx)%equationsSetIndex=equationsSetIdx
                solverMapping%equationsSetToSolverMap(equationsSetIdx)%solverMapping=>solverMapping
                solverMapping%equationsSetToSolverMap(equationsSetIdx)%equations=>equations
                dependentField=>equationsSet%dependent%DEPENDENT_FIELD
                IF(ASSOCIATED(dependentField)) THEN
                  lhsVariableType=solverMapping%createValuesCache%lhsVariableType(equationsSetIdx)
                  lhsVariable=>dependentField%VARIABLE_TYPE_MAP(lhsVariableType)%ptr
                  IF(ASSOCIATED(lhsVariable)) THEN
                    !See if this variable is already in the list of row variables.
                    CALL SolverMapping_SolverVariableFindAdd(solverMapping%rowVariablesList,lhsVariable, &
                      & SOLVER_MAPPING_EQUATIONS_EQUATIONS_SET,equationsSetIdx,err,error,*999)                    
                    !Find column list variables
                    equationsMapping=>equations%EQUATIONS_MAPPING
                    IF(ASSOCIATED(equationsMapping)) THEN
                      dynamicMapping=>equationsMapping%DYNAMIC_MAPPING
                      linearMapping=>equationsMapping%LINEAR_MAPPING
                      nonlinearMapping=>equationsMapping%NONLINEAR_MAPPING
                      DO solverMatrixIdx=1,solverMapping%numberOfSolverMatrices
                        IF(ASSOCIATED(dynamicMapping)) THEN
                          dynamicVariableType=solverMapping%createValuesCache%dynamicVariableType(equationsSetIdx)
                          dynamicVariable=>dependentField%VARIABLE_TYPE_MAP(dynamicVariableType)%ptr
                          IF(ASSOCIATED(dynamicVariable)) THEN
                            !See if this variable is already in the list of column variables.
                            CALL SolverMapping_SolverVariableFindAdd(solverMapping%columnVariablesList(solverMatrixIdx), &
                              & dynamicVariable,SOLVER_MAPPING_EQUATIONS_EQUATIONS_SET,equationsSetIdx,err,error,*999)
                          ELSE
                            localError="Dependent field variable is not associated for dynamic variable type "// &
                              & TRIM(NumberToVstring(dynamicVariableType,"*",err,error))//" in equations set index "// &
                              & TRIM(NumberToVstring(equationsSetIdx,"*",err,error))//"."
                            CALL FlagError(localError,err,error,*999)
                          ENDIF
                        ENDIF
                        IF(ASSOCIATED(linearMapping)) THEN
                          DO variableIdx=1,solverMapping%createValuesCache%matrixVariableTypes(0,equationsSetIdx,solverMatrixIdx)
                            linearVariableType=solverMapping%createValuesCache%matrixVariableTypes(variableIdx,equationsSetIdx, &
                              & solverMatrixIdx)
                            linearVariable=>dependentField%VARIABLE_TYPE_MAP(linearVariableType)%ptr
                            IF(ASSOCIATED(linearVariable)) THEN
                              !See if this variable is already in the list of column variables.
                              CALL SolverMapping_SolverVariableFindAdd(solverMapping%columnVariablesList(solverMatrixIdx), &
                                & linearVariable,SOLVER_MAPPING_EQUATIONS_EQUATIONS_SET,equationsSetIdx,err,error,*999)
                            ELSE
                              localError="Dependent field variable is not associated for linear variable type "// &
                                & TRIM(NumberToVstring(linearVariableType,"*",err,error))//" in equations set index "// &
                                & TRIM(NumberToVstring(equationsSetIdx,"*",err,error))//"."
                              CALL FlagError(localError,err,error,*999)
                            ENDIF
                          ENDDO !variableIdx
                        ENDIF
                        IF(ASSOCIATED(nonlinearMapping)) THEN
                          DO variableIdx=1,solverMapping%createValuesCache%residualVariableTypes(0,equationsSetIdx,solverMatrixIdx)
                            residualVariableType=solverMapping%createValuesCache%matrixVariableTypes(variableIdx,equationsSetIdx, &
                              & solverMatrixIdx)
                            residualVariable=>dependentField%VARIABLE_TYPE_MAP(residualVariableType)%ptr
                            IF(ASSOCIATED(residualVariable)) THEN
                              !See if this variable is already in the list of column variables.
                              CALL SolverMapping_SolverVariableFindAdd(solverMapping%columnVariablesList(solverMatrixIdx), &
                                & residualVariable,SOLVER_MAPPING_EQUATIONS_EQUATIONS_SET,equationsSetIdx,err,error,*999)
                            ELSE
                              localError="Dependent field variable is not associated for residual variable type "// &
                                & TRIM(NumberToVstring(residualVariableType,"*",err,error))//" in equations set index "// &
                                & TRIM(NumberToVstring(equationsSetIdx,"*",err,error))//"."
                              CALL FlagError(localError,err,error,*999)
                            ENDIF
                          ENDDO !variableIdx
                        ENDIF
                      ENDDO !solverMatrixIdx
                    ELSE
                      CALL FlagError("Equations mapping is not associated.",err,error,*999)
                    ENDIF
                  ELSE
                    localError="Dependent field variable is not associated for LHS variable type "// &
                      & TRIM(NumberToVstring(lhsVariableType,"*",err,error))//" in equations set index "// &
                      & TRIM(NumberToVstring(equationsSetIdx,"*",err,error))//"."
                    CALL FlagError(localError,err,error,*999)
                  ENDIF
                ELSE
                  CALL FlagError("Equations set dependent field is not associated.",err,error,*999)
                ENDIF
              ELSE
                CALL FlagError("Equations set equations is not associated.",err,error,*999)
              ENDIF
            ELSE
              CALL FlagError("Equations set is not associated.",err,error,*999)
            ENDIF
          ENDDO !equationsSetIdx
          DO interfaceConditionIdx=1,solverMapping%numberOfInterfaceConditions
            interfaceCondition=>solverMapping%interfaceConditions(interfaceConditionIdx)%PTR
            IF(ASSOCIATED(interfaceCondition)) THEN
              interfaceEquations=>interfaceCondition%INTERFACE_EQUATIONS
              IF(ASSOCIATED(interfaceEquations)) THEN
                CALL SolverMapping_InterfaceConditionToSolverMapInitialise(solverMapping%interfaceConditionToSolverMap( &
                  & interfaceConditionIdx),err,error,*999)
                solverMapping%interfaceConditionToSolverMap(interfaceConditionIdx)%interfaceConditionIndex= &
                  & interfaceConditionIdx
                solverMapping%interfaceConditionToSolverMap(interfaceConditionIdx)%solverMapping=>solverMapping
                solverMapping%interfaceConditionToSolverMap(interfaceConditionIdx)%interfaceEquations=>interfaceEquations
              ELSE
                CALL FlagError("Interface condition interface equations is not associated.",err,error,*999)
              ENDIF
            ELSE
              CALL FlagError("Interface condition is not associated.",err,error,*999)
            ENDIF
          ENDDO !interfaceConditionIdx
                 
          !
          ! Allocate and initialise
          !
          !Compute the order of the row variables for the solver matrices          
          CALL List_DetachAndDestroy(solverMapping%createValuesCache%equationsRowVariablesList%ptr, &
            & numberOfEquationsRowVariables,equationsRowVariables,err,error,*999)
          CALL List_DetachAndDestroy(solverMapping%createValuesCache%interfaceRowVariablesList%ptr, &
            & numberOfInterfaceRowVariables,interfaceRowVariables,err,error,*999)

          variableProcessed=.FALSE.
          solverVariableIdx=0
          DO variableIdx=1,numberOfEquationsRowVariables

          DO equationsSetIdx=1,solverMapping%numberOfEquationsSets
            equationsSet=>solverMapping%equationsSets(equationsSetIdx)%ptr
            IF(ASSOCIATED(equationsSet)) THEN
              equations=>equationsSet%equations
              IF(ASSOCIATED(equations)) THEN
                 CALL SolverMapping_EquationsSetToSolverMapInitialise(solverMapping%equationsSetToSolverMap( &
                  & equationsSetIdx),err,error,*999)
                solverMapping%equationsSetToSolverMap(equationsSetIdx)%equationsSetIndex=equationsSetIdx
                solverMapping%equationsSetToSolverMap(equationsSetIdx)%solverMapping=>solverMapping
                solverMapping%equationsSetToSolverMap(equationsSetIdx)%equations=>equations
                
                
           !Calculate the column mappings for each solver matrix
          DO solverMatrixIdx=1,solverMapping%numberOfSolverMatrices

            !Initialise the column variables list
            CALL SolverMapping_VariablesInitialise(solverMapping%columnVariablesList(solverMatrixIdx),err,error,*999)
            !
            ! 4a Calculate the list of field variables involved in the columns of the solver matrix
            !
            !Compute the order of variables for the solver matrices          
            CALL List_DetachAndDestroy(solverMapping%createValuesCache%equationsColVariablesList(solverMatrixIdx)%ptr, &
              & numberOfEquationsVariables,equationsVariables,err,error,*999)
            CALL List_DetachAndDestroy(solverMapping%createValuesCache%interfaceColVariablesList(solverMatrixIdx)%ptr, &
              & numberOfInterfaceVariables,interfaceVariables,err,error,*999)
            ALLOCATE(solverMapping%columnVariablesList(solverMatrixIdx)%variables(numberOfEquationsVariables+ &
              & numberOfInterfaceVariables),STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate column variables list variables.",err,error,*999)
            solverMapping%columnVariablesList(solverMatrixIdx)%numberOfVariables=numberOfEquationsVariables+ &
              & numberOfInterfaceVariables
            ALLOCATE(variablesList(numberOfEquationsVariables+numberOfInterfaceVariables),STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate variables list.",err,error,*999)
            ALLOCATE(variableProcessed(numberOfEquationsVariables+numberOfInterfaceVariables),STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate variable processed.",err,error,*999)
            variableProcessed=.FALSE.
            solverVariableIdx=0
            DO variableIdx=1,numberOfEquationsVariables
              solverVariableIdx=solverVariableIdx+1
              CALL SolverMapping_VariableInitialise(solverMapping%columnVariablesList(solverMatrixIdx)% &
                & variables(solverVariableIdx),err,error,*999)
              equationsSetIdx=equationsVariables(1,variableIdx)
              equationsSet=>solverMapping%equationsSets(equationsSetIdx)%ptr
              IF(ASSOCIATED(equationsSet)) THEN
                dependentField=>equationsSet%dependent%DEPENDENT_FIELD
                IF(ASSOCIATED(dependentField)) THEN
                  variableType=equationsVariables(2,variableIdx)
                  variable=>dependentField%VARIABLE_TYPE_MAP(variableType)%ptr
                  IF(ASSOCIATED(variable)) THEN
                    solverMapping%columnVariablesList(solverMatrixIdx)%variables(solverVariableIdx)%variable=>variable
                    solverMapping%columnVariablesList(solverMatrixIdx)%variables(solverVariableIdx)%variableType=variableType
                    NULLIFY(variablesList(solverVariableIdx)%ptr)
                    CALL List_CreateStart(variablesList(solverVariableIdx)%ptr,err,error,*999)
                    CALL List_DataTypeSet(variablesList(solverVariableIdx)%ptr,LIST_INTG_TYPE,err,error,*999)
                    CALL List_DataDimensionSet(variablesList(solverVariableIdx)%ptr,3,err,error,*999)
                    CALL List_KeyDimensionSet(variablesList(solverVariableIdx)%ptr,1,err,error,*999)
                    CALL List_CreateFinish(variablesList(solverVariableIdx)%ptr,err,error,*999)
                  ELSE
                    CALL FlagError("Dependent field variable is not associated.",err,error,*999)
                  ENDIF
                ELSE
                  CALL FlagError("Equations set dependent field is not associated.",err,error,*999)
                ENDIF
              ELSE
                CALL FlagError("Equations set is not associated.",err,error,*999)
              ENDIF
            ENDDO !variableIdx
            IF(ALLOCATED(equationsVariables)) DEALLOCATE(equationsVariables)
            DO variableIdx=1,numberOfInterfaceVariables
              solverVariableIdx=solverVariableIdx+1
              CALL SolverMapping_VariableInitialise(solverMapping%columnVariablesList(solverMatrixIdx)% &
                & variables(solverVariableIdx),err,error,*999)
              interfaceConditionIdx=interfaceVariables(1,variableIdx)
              interfaceCondition=>solverMapping%interfaceConditions(interfaceConditionIdx)%ptr
              IF(ASSOCIATED(interfaceCondition)) THEN
                IF(ASSOCIATED(interfaceCondition%lagrange)) THEN
                  lagrangeField=>interfaceCondition%lagrange%LAGRANGE_FIELD
                  IF(ASSOCIATED(lagrangeField)) THEN
                    variableType=interfaceVariables(2,variableIdx)
                    variable=>lagrangeField%VARIABLE_TYPE_MAP(variableType)%ptr
                    IF(ASSOCIATED(variable)) THEN
                      solverMapping%columnVariablesList(solverMatrixIdx)%variables(solverVariableIdx)%variable=>variable
                      solverMapping%columnVariablesList(solverMatrixIdx)%variables(solverVariableIdx)%variableType=variableType
                      NULLIFY(variablesList(solverVariableIdx)%ptr)
                      CALL List_CreateStart(variablesList(solverVariableIdx)%ptr,err,error,*999)
                      CALL List_DataTypeSet(variablesList(solverVariableIdx)%ptr,LIST_INTG_TYPE,err,error,*999)
                      CALL List_DataDimensionSet(variablesList(solverVariableIdx)%ptr,3,err,error,*999)
                      CALL List_KeyDimensionSet(variablesList(solverVariableIdx)%ptr,1,err,error,*999)
                      CALL List_CreateFinish(variablesList(solverVariableIdx)%ptr,err,error,*999)
                    ELSE
                      CALL FlagError("Lagrange field variable is not associated.",err,error,*999)
                    ENDIF
                  ELSE
                    CALL FlagError("Interface condition Lagrange field is not associated.",err,error,*999)
                  ENDIF
                ELSE
                  CALL FlagError("Interface condition Lagrange is not associated.",err,error,*999)
                ENDIF
              ELSE
                CALL FlagError("Interface condition is not associated.",err,error,*999)
              ENDIF
            ENDDO !variableIdx
            IF(ALLOCATED(interfaceVariables)) DEALLOCATE(interfaceVariables)
           
            !Allocate sub-matrix information
            !SUB_MATRIX_INFORMATION(1,equations_idx,variable_idx) = The equations type, see SOLVER_MAPPING_EquationsTypes
            !SUB_MATRIX_INFORMATION(2,equations_idx,variable_idx) = The equations set or interface condition index
            !SUB_MATRIX_INFORMATION(3,equations_idx,variable_idx) = The interface matrix index, or 0 for an equations set matrix
            !equations_idx goes from 1 to the number of equations sets + interface conditions
            !variable_idx goes from 1 to the number of variables mapped to this solver matrix
             ALLOCATE(subMatrixInformation(3,solverMapping%numberOfEquationsSets+solverMapping%numberOfInterfaceConditions, &
              & solverMapping%columnVariablesList(solverMatrixIdx)%numberOfVariables),STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate sub matrix information.",err,error,*999)
            subMatrixInformation=0
            !Allocate sub-matrix list information
            ALLOCATE(subMatrixList(0:3,solverMapping%numberOfEquationsSets+solverMapping%numberOfInterfaceConditions, &
              & solverMapping%columnVariablesList(solverMatrixIdx)%numberOfVariables),STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate sub matrix list.",err,error,*999)
            subMatrixList=0

          
          ALLOCATE(interfaceEquationsLists(solverMapping%numberOfInterfaceConditions),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate equations set list.",err,error,*999)
          DO interfaceConditionIdx=1,solverMapping%numberOfInterfaceConditions
            interfaceCondition=>solverMapping%interfaceConditions(interfaceConditionIdx)%PTR
            IF(ASSOCIATED(interfaceCondition)) THEN
              interfaceEquations=>interfaceCondition%INTERFACE_EQUATIONS
              IF(ASSOCIATED(interfaceEquations)) THEN
                CALL SolverMapping_InterfaceConditionToSolverMapInitialise(solverMapping%interfaceConditionToSolverMap( &
                  & interfaceConditionIdx),err,error,*999)
                solverMapping%interfaceConditionToSolverMap(interfaceConditionIdx)%interfaceConditionIndex= &
                  & interfaceConditionIdx
                solverMapping%interfaceConditionToSolverMap(interfaceConditionIdx)%solverMapping=>solverMapping
                solverMapping%interfaceConditionToSolverMap(interfaceConditionIdx)%interfaceEquations=>interfaceEquations
                NULLIFY(interfaceEquationsLists(interfaceConditionIdx)%ptr)
                CALL List_CreateStart(interfaceEquationsLists(interfaceConditionIdx)%ptr,err,error,*999)
                CALL List_DataTypeSet(interfaceEquationsLists(interfaceConditionIdx)%ptr,LIST_INTG_TYPE,err,error,*999)
                CALL List_DataDimensionSet(interfaceEquationsLists(interfaceConditionIdx)%ptr,2,err,error,*999)
                CALL List_CreateFinish(interfaceEquationsLists(interfaceConditionIdx)%ptr,err,error,*999)
              ELSE
                CALL FlagError("Interface condition interface equations is not associated.",err,error,*999)
              ENDIF
            ELSE
              CALL FlagError("Interface condition is not associated.",err,error,*999)
            ENDIF
          ENDDO !interfaceConditionIdx
          !
          ! Loop over equations sets
          !
          DO equationsSetIdx=1,solverMapping%numberOfEquationsSets
            equationsSet=>solverMapping%equationsSets(equationsSetIdx)%ptr
            IF(ASSOCIATED(equationsSet)) THEN
              equations=>equationsSet%equations
              IF(ASSOCIATED(equations)) THEN
                CALL SolverMapping_EquationsSetToSolverMapInitialise(solverMapping%equationsSetToSolverMap( &
                  & equationsSetIdx),err,error,*999)
                solverMapping%equationsSetToSolverMap(equationsSetIdx)%equationsSetIndex=equationsSetIdx
                solverMapping%equationsSetToSolverMap(equationsSetIdx)%solverMapping=>solverMapping
                solverMapping%equationsSetToSolverMap(equationsSetIdx)%equations=>equations
                !Set up list of interface conditions affecting this equations set
                CALL List_DetachAndDestroy(solverMapping%createValuesCache%interfaceIndices(equationsSetIdx)%ptr, &
                  & numberOfInterfaces,interfaceEquationsList,err,error,*999)
                ALLOCATE(solverMapping%equationsSetToSolverMap(equationsSetIdx)% &
                  & equationsToSolverMatrixMapsInterface(numberOfInterfaces),STAT=err)
                IF(err/=0) CALL FlagError("Could not allocate equations to solver maps interface.",err,error,*999)
                solverMapping%equationsSetToSolverMap(equationsSetIdx)%numberOfInterfaceConditions=numberOfInterfaces
                DO interfaceIdx=1,numberOfInterfaces
                  CALL SolverMapping_EquationsToSolverInterfaceInitialise(solverMapping%equationsSetToSolverMap( &
                    & equationsSetIdx)%equationsToSolverMatrixMapsInterface(interfaceIdx),err,error,*999)
                  interfaceConditionIdx=interfaceEquationsList(1,interfaceIdx)
                  interfaceMatrixIdx=interfaceEquationsList(2,interfaceIdx)
                  interfaceCondition=>solverMapping%interfaceConditions(interfaceConditionIdx)%ptr
                  IF(ASSOCIATED(interfaceCondition)) THEN
                    solverMapping%equationsSetToSolverMap(equationsSetIdx)%equationsToSolverMatrixMapsInterface( &
                      & interfaceIdx)%interfaceConditionIndex=interfaceConditionIdx
                    solverMapping%equationsSetToSolverMap(equationsSetIdx)%equationsToSolverMatrixMapsInterface( &
                      & interfaceIdx)%interfaceCondition=>interfaceCondition
                    solverMapping%equationsSetToSolverMap(equationsSetIdx)%equationsToSolverMatrixMapsInterface( &
                      & interfaceIdx)%interfaceMatrixNumber=interfaceMatrixIdx
                    interfaceEquationsListItem(1)=equationsSetIdx
                    interfaceEquationsListItem(2)=interfaceMatrixIdx
                    CALL List_ItemAdd(interfaceEquationsLists(interfaceConditionIdx)%ptr,interfaceEquationsListItem,err,error,*999)
                  ELSE
                    CALL FlagError("Interface condition is not associated.",err,error,*999)
                  ENDIF
                ENDDO !interfaceConditionIdx
                IF(ALLOCATED(interfaceEquationsList)) DEALLOCATE(interfaceEquationsList)
              ELSE
                CALL FlagError("Equations set equations is not associated.",err,error,*999)
              ENDIF
            ELSE
              CALL FlagError("Equations set is not associated.",err,error,*999)
            ENDIF
          ENDDO !equationsSetIdx
          DO interfaceConditionIdx=1,solverMapping%numberOfInterfaceConditions
            CALL List_DetachAndDestroy(interfaceEquationsLists(interfaceConditionIdx)%ptr,numberOfEquationsSets, &
              & interfaceEquationsList,err,error,*999)
            ALLOCATE(solverMapping%interfaceConditionToSolverMap(interfaceConditionIdx)% &
              & interfaceToSolverMatrixMapsEquations(numberOfEquationsSets),STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate interface to solver maps equations.",err,error,*999)
            solverMapping%interfaceConditionToSolverMap(interfaceConditionIdx)%numberOfEquationsSets=numberOfEquationsSets
            DO equationsIdx=1,numberOfEquationsSets
              CALL SolverMapping_InterfaceToSolverEquationsInitialise(solverMapping%interfaceConditionToSolverMap( &
                & interfaceConditionIdx)%interfaceToSolverMatrixMapsEquations(equationsIdx),err,error,*999)
              equationsSetIdx=interfaceEquationsList(1,equationsIdx)
              interfaceMatrixIdx=interfaceEquationsList(2,equationsIdx)
              equationsSet=>solverMapping%equationsSets(equationsSetIdx)%ptr
              IF(ASSOCIATED(equationsSet)) THEN
                solverMapping%interfaceConditionToSolverMap(interfaceConditionIdx)% &
                  & interfaceToSolverMatrixMapsEquations(equationsIdx)%equationsSetIndex=equationsSetIdx
                solverMapping%interfaceConditionToSolverMap(interfaceConditionIdx)% &
                  & interfaceToSolverMatrixMapsEquations(equationsIdx)%equationsSet=>equationsSet
                solverMapping%interfaceConditionToSolverMap(interfaceConditionIdx)% &
                  & interfaceToSolverMatrixMapsEquations(equationsIdx)%interfaceMatrixIndex=interfaceMatrixIdx
              ELSE
                CALL FlagError("Equations set is not associated.",err,error,*999)
              ENDIF
            ENDDO !equationsIdx
            IF(ALLOCATED(interfaceEquationsList)) DEALLOCATE(interfaceEquationsList)
          ENDDO !interfaceConditionIdx
          !
          !--- Row mappings ---
          !
          ! 2. Determine the number of rows in the solver matrix. Do this the by setting up a list of rows for each rank.
          !    We can then later arrange the rows in rank order by looping over the ranks in the list and then the rows
          !    for each rank.
          !
          !Calculate the row mappings.
          !We do not have any couplings defined at the moment there is only a 1-1 mapping.
          myRank=COMPUTATIONAL_ENVIRONMENT%MY_COMPUTATIONAL_NODE_NUMBER
          numberOfGlobalSolverRows=0
          numberOfLocalSolverRows=0
          !Add in the rows from any equations sets that have been added to the solver equations
          !Presort the row numbers by rank.
          !
          !Initialise the row variables list
          CALL SolverMapping_VariablesInitialise(solverMapping%rowVariablesList,err,error,*999)
          !
          ! 2a Calculate the list of field variables involved in the rows of the solver matrix
          !
          !Compute the order of variables      
          CALL List_DetachAndDestroy(solverMapping%createValuesCache%equationsRowVariablesList,numberOfEquationsVariables, &
            & equationsVariables,err,error,*999)
          CALL List_DetachAndDestroy(solverMapping%createValuesCache%interfaceRowVariablesList,numberOfInterfaceVariables, &
            & interfaceVariables,err,error,*999)
          ALLOCATE(solverMapping%rowVariablesList%variables(numberOfEquationsVariables+numberOfInterfaceVariables),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate row variables list variables.",err,error,*999)          
          solverMapping%rowVariablesList%numberOfVariables=numberOfEquationsVariables+numberOfInterfaceVariables
          ALLOCATE(variablesList(numberOfEquationsVariables+numberOfInterfaceVariables),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate variables list.",err,error,*999)
          ALLOCATE(variableProcessed(numberOfEquationsVariables+numberOfInterfaceVariables),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate variable processed.",err,error,*999)
          variableProcessed=.FALSE.
          solverVariableIdx=0
          DO variableIdx=1,numberOfEquationsVariables
            solverVariableIdx=solverVariableIdx+1
            CALL SolverMapping_VariableInitialise(solverMapping%rowVariablesList%variables(solverVariableIdx),err,error,*999)
            equationsSetIdx=equationsVariables(1,variableIdx)
            equationsSet=>solverMapping%equationsSets(equationsSetIdx)%ptr
            IF(ASSOCIATED(equationsSet)) THEN
              dependentField=>equationsSet%dependent%DEPENDENT_FIELD
              IF(ASSOCIATED(dependentField)) THEN
                variableType=equationsVariables(2,variableIdx)
                variable=>dependentField%VARIABLE_TYPE_MAP(variableType)%ptr
                IF(ASSOCIATED(variable)) THEN
                  solverMapping%rowVariablesList%variables(solverVariableIdx)%variable=>variable
                  solverMapping%rowVariablesList%variables(solverVariableIdx)%variableType=variableType
                  NULLIFY(variablesList(solverVariableIdx)%ptr)
                  CALL List_CreateStart(variablesList(solverVariableIdx)%ptr,err,error,*999)
                  CALL List_DataTypeSet(variablesList(solverVariableIdx)%ptr,LIST_INTG_TYPE,err,error,*999)
                  CALL List_DataDimensionSet(variablesList(solverVariableIdx)%ptr,3,err,error,*999)
                  CALL List_KeyDimensionSet(variablesList(solverVariableIdx)%ptr,1,err,error,*999)
                  CALL List_CreateFinish(variablesList(solverVariableIdx)%ptr,err,error,*999)
                ELSE
                  CALL FlagError("Dependent field variable is not associated.",err,error,*999)
                ENDIF
              ELSE
                CALL FlagError("Equations set dependent field is not associated.",err,error,*999)
              ENDIF
            ELSE
              CALL FlagError("Equations set is not associated.",err,error,*999)
            ENDIF
          ENDDO !variableIdx
          IF(ALLOCATED(equationsVariables)) DEALLOCATE(equationsVariables)
          DO variableIdx=1,numberOfInterfaceVariables
            solverVariableIdx=solverVariableIdx+1
            CALL SolverMapping_VariableInitialise(solverMapping%rowVariablesList%variables(solverVariableIdx),err,error,*999)
            interfaceConditionIdx=interfaceVariables(1,variableIdx)
            interfaceCondition=>solverMapping%interfaceConditions(interfaceConditionIdx)%ptr
            IF(ASSOCIATED(interfaceCondition)) THEN
              IF(ASSOCIATED(interfaceCondition%lagrange)) THEN
                lagrangeField=>interfaceCondition%lagrange%LAGRANGE_FIELD
                IF(ASSOCIATED(lagrangeField)) THEN
                  variableType=interfaceVariables(2,variableIdx)
                  variable=>lagrangeField%VARIABLE_TYPE_MAP(variableType)%ptr
                  IF(ASSOCIATED(variable)) THEN
                    solverMapping%columnVariablesList(solverMatrixIdx)%variables(solverVariableIdx)%variable=>variable
                    solverMapping%columnVariablesList(solverMatrixIdx)%variables(solverVariableIdx)%variableType=variableType
                    NULLIFY(variablesList(solverVariableIdx)%ptr)
                    CALL List_CreateStart(variablesList(solverVariableIdx)%ptr,err,error,*999)
                    CALL List_DataTypeSet(variablesList(solverVariableIdx)%ptr,LIST_INTG_TYPE,err,error,*999)
                    CALL List_DataDimensionSet(variablesList(solverVariableIdx)%ptr,3,err,error,*999)
                    CALL List_KeyDimensionSet(variablesList(solverVariableIdx)%ptr,1,err,error,*999)
                    CALL List_CreateFinish(variablesList(solverVariableIdx)%ptr,err,error,*999)
                  ELSE
                    CALL FlagError("Lagrange field variable is not associated.",err,error,*999)
                  ENDIF
                ELSE
                  CALL FlagError("Interface condition Lagrange field is not associated.",err,error,*999)
                ENDIF
              ELSE
                CALL FlagError("Interface condition Lagrange is not associated.",err,error,*999)
              ENDIF
            ELSE
              CALL FlagError("Interface condition is not associated.",err,error,*999)
            ENDIF
          ENDDO !variableIdx
          IF(ALLOCATED(interfaceVariables)) DEALLOCATE(interfaceVariables)

          !Allocate and initialise the rank lists.
          ALLOCATE(rankGlobalRowsLists(solverMapping%numberOfEquationsSets+solverMapping%numberOfInterfaceConditions, &
            & 0:COMPUTATIONAL_ENVIRONMENT%NUMBER_COMPUTATIONAL_NODES-1),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate rank global rows lists.",err,error,*999)
          CALL SolverDofCouplings_Initialise(rowCouplings,err,error,*999)
          DO rank=0,COMPUTATIONAL_ENVIRONMENT%NUMBER_COMPUTATIONAL_NODES-1
            equationsIdx=0
            DO equationsSetIdx=1,solverMapping%numberOfEquationsSets
              equationsIdx=equationsIdx+1
              equationsSet=>solverMapping%equationsSets(equationsSetIdx)%ptr
              IF(ASSOCIATED(equationsSet)) THEN
                equations=>equationsSet%equations
                IF(ASSOCIATED(equations)) THEN
                  equationsMapping=>equations%EQUATIONS_MAPPING                
                  IF(ASSOCIATED(equationsMapping)) THEN
                    NULLIFY(rankGlobalRowsLists(equationsIdx,rank)%ptr)
                    CALL List_CreateStart(rankGlobalRowsLists(equationsIdx,rank)%ptr,err,error,*999)
                    CALL List_DataTypeSet(rankGlobalRowsLists(equationsIdx,rank)%ptr,LIST_INTG_TYPE,err,error,*999)
                    CALL List_InitialSizeSet(rankGlobalRowsLists(equationsIdx,rank)%ptr,INT(equationsMapping% &
                      & NUMBER_OF_GLOBAL_ROWS/COMPUTATIONAL_ENVIRONMENT%NUMBER_COMPUTATIONAL_NODES,INTG), &
                      & err,error,*999)
                    CALL List_DataDimensionSet(rankGlobalRowsLists(equationsIdx,rank)%ptr,4,err,error,*999)
                    CALL List_KeyDimensionSet(rankGlobalRowsLists(equationsIdx,rank)%ptr,1,err,error,*999)
                    CALL List_CreateFinish(rankGlobalRowsLists(equationsIdx,rank)%ptr,err,error,*999)                    
                  ELSE
                    CALL FlagError("Equations equations mapping is not associated",err,error,*999)
                  ENDIF
                ELSE
                  CALL FlagError("Equations set equations is not associated.",err,error,*999)
                ENDIF
              ELSE
                CALL FlagError("Equations set is not associated.",err,error,*999)
              ENDIF
            ENDDO !equationsSetIdx
            DO interfaceConditionIdx=1,solverMapping%numberOfInterfaceConditions
              equationsIdx=equationsIdx+1
              interfaceCondition=>solverMapping%interfaceConditions(interfaceConditionIdx)%ptr
              IF(ASSOCIATED(interfaceCondition)) THEN
                SELECT CASE(interfaceCondition%method)
                CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD)
                 
                  interfaceEquations=>interfaceCondition%INTERFACE_EQUATIONS
                  IF(ASSOCIATED(interfaceEquations)) THEN
                    interfaceMapping=>interfaceEquations%INTERFACE_MAPPING
                    IF(ASSOCIATED(interfaceMapping)) THEN
                      NULLIFY(rankGlobalRowsLists(equationsIdx,rank)%ptr)
                      CALL List_CreateStart(rankGlobalRowsLists(equationsIdx,rank)%ptr,err,error,*999)
                      CALL List_DataTypeSet(rankGlobalRowsLists(equationsIdx,rank)%ptr,LIST_INTG_TYPE,err,error,*999)
                      CALL List_InitialSizeSet(rankGlobalRowsLists(equationsIdx,rank)%ptr, &
                        & INT(interfaceMapping%NUMBER_OF_GLOBAL_COLUMNS/COMPUTATIONAL_ENVIRONMENT%NUMBER_COMPUTATIONAL_NODES, &
                        & INTG),err,error,*999)
                      CALL List_DataDimensionSet(rankGlobalRowsLists(equationsIdx,rank)%ptr,4,err,error,*999)
                      CALL List_KeyDimensionSet(rankGlobalRowsLists(equationsIdx,rank)%ptr,1,err,error,*999)
                      CALL List_CreateFinish(rankGlobalRowsLists(equationsIdx,rank)%ptr,err,error,*999)                      
                    ELSE
                      CALL FlagError("Interface equations interface mapping is not associated.",err,error,*999)
                    ENDIF
                  ELSE
                    CALL FlagError("Interface condition interface equations is not associated.",err,error,*999)
                  ENDIF
                CASE(INTERFACE_CONDITION_AUGMENTED_LAGRANGE_METHOD)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(INTERFACE_CONDITION_PENALTY_METHOD)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(INTERFACE_CONDITION_POINT_TO_POINT_METHOD)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE DEFAULT
                  localError="The interface condition method of "// &
                    & TRIM(NumberToVString(interfaceCondition%method,"*",err,error))//" is invalid."
                  CALL FlagError(localError,err,error,*999)
                END SELECT                
              ELSE
                CALL FlagError("Interface condition is not associated.",err,error,*999)
              ENDIF
            ENDDO !interfaceConditionIdx
          ENDDO !rank

          !
          ! 2b Calculate the number of rows
          !          
          
          !Calculate the number of local and global rows. Do this by looking at the boundary conditions for the LHS field variable
          !in the row. If the variable is set as a fixed boundary condition then do not include the row.
          
          !Calculate the number of solver dofs
          ALLOCATE(numberOfVariableGlobalSolverDofs(solverMapping%rowVariablesList%numberOfVariables),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate number of global solver dofs.",err,error,*999)
          ALLOCATE(numberOfVariableLocalSolverDofs(solverMapping%rowVariablesList%numberOfVariables),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate number of local solver dofs.",err,error,*999)
          ALLOCATE(totalNumberOfVariableLocalSolverDofs(solverMapping%rowVariablesList%numberOfVariables),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate total number of local solver dofs.",err,error,*999)
          
          numberOfVariableGlobalSolverDofs=0
          numberOfVariableLocalSolverDofs=0
          totalNumberOfVariableLocalSolverDofs=0
            
          equationsIdx=0
          DO equationsSetIdx=1,solverMapping%numberOfEquationsSets
            equationsIdx=equationsIdx+1
            equationsSet=>solverMapping%equationsSets(equationsSetIdx)%ptr
            IF(ASSOCIATED(equationsSet)) THEN
              equations=>equationsSet%equations
              IF(ASSOCIATED(equations)) THEN
                equationsMapping=>equations%EQUATIONS_MAPPING                
                IF(ASSOCIATED(equationsMapping)) THEN
                  lhsMapping=>equationsMapping%lhsMapping
                  IF(ASSOCIATED(lhsMapping)) THEN
                    rowDofsMapping=>equationsMapping%ROW_DOFS_MAPPING
                    IF(ASSOCIATED(rowDofsMapping)) THEN                    
                      dependentField=>equationsSet%dependent%DEPENDENT_FIELD
                      IF(ASSOCIATED(dependentField)) THEN
                        boundaryConditions=>solverEquations%BOUNDARY_CONDITIONS
                        IF(ASSOCIATED(boundaryConditions)) THEN
                          !Loop over the global rows for this equations set
                          DO globalRow=1,equationsMapping%NUMBER_OF_GLOBAL_ROWS
                            !Find the rank that owns this global row
                            rowRank=-1
                            DO rankIdx=1,rowDofsMapping%GLOBAL_TO_LOCAL_MAP(globalRow)%NUMBER_OF_DOMAINS
                              IF(rowDofsMapping%GLOBAL_TO_LOCAL_MAP(globalRow)%LOCAL_TYPE(rankIdx)/=DOMAIN_LOCAL_GHOST) THEN
                                rowRank=rowDofsMapping%GLOBAL_TO_LOCAL_MAP(globalRow)%DOMAIN_NUMBER(rankIdx)
                                localRow=rowDofsMapping%GLOBAL_TO_LOCAL_MAP(globalRow)%LOCAL_NUMBER(rankIdx)
                                EXIT
                              ENDIF
                            ENDDO !rankIdx
                            IF(rowRank>=0) THEN
                              includeRow=.TRUE.
                              constrainedDof=.FALSE.
                              globalDofCouplingNumber=0
                              dependentVariable=>lhsMapping%lhsVariable
                              CALL BOUNDARY_CONDITIONS_VARIABLE_GET(boundaryConditions,dependentVariable, &
                                & boundaryConditionsVariable,err,error,*999)
                              IF(ASSOCIATED(boundaryConditionsVariable)) THEN
                                !This is wrong as we only have the mappings for the local rank not the global ranks.
                                !For now assume 1-1 mapping between rows and dofs.
                                globalDof=globalRow
                                includeRow=includeRow.AND.(boundaryConditionsVariable%DOF_TYPES(globalDof)== &
                                  & BOUNDARY_CONDITION_DOF_FREE)
                                constrainedDof=constrainedDof.OR.(boundaryConditionsVariable%DOF_TYPES(globalDof)== &
                                  & BOUNDARY_CONDITION_DOF_CONSTRAINED)
                                dofConstraints=>boundaryConditionsVariable%dofConstraints
                                IF(ASSOCIATED(dofConstraints)) THEN
                                  IF(dofConstraints%numberOfConstraints>0) THEN
                                    IF(ALLOCATED(dofConstraints%dofCouplings)) THEN
                                      IF(ASSOCIATED(dofConstraints%dofCouplings(globalDof)%ptr)) THEN
                                        !This equations row is the owner of a solver row that is mapped to
                                        !multiple other equations rows, add it to the list of global row
                                        !couplings and remember the index into the global list for this solver row
                                        CALL SolverDofCouplings_AddCoupling(rowCouplings, &
                                          & dofConstraints%dofCouplings(globalDof)%ptr, &
                                          & globalDofCouplingNumber,err,error,*999)
                                      END IF
                                    ELSE
                                      CALL FlagError("DOF constraints DOF couplings are not allocated.",err,error,*999)
                                    END IF
                                  END IF
                                END IF
                              ELSE
                                CALL FlagError("Boundary condition variable is not associated.",err,error,*999)
                              ENDIF
                              rowListItem(1)=globalRow
                              rowListItem(2)=localRow
                              IF(includeRow) THEN
                                rowListItem(3)=1
                                numberOfGlobalSolverRows=numberOfGlobalSolverRows+1
                                !Don't need to worry about ghosted rows.
                                IF(rowRank==myrank) numberOfLocalSolverRows=numberOfLocalSolverRows+1 !1-1 mapping
                              ELSE IF(constrainedDof) THEN
                                rowListItem(3)=2
                              ELSE
                                rowListItem(3)=0
                              ENDIF !include row
                              rowListItem(4)=globalDofCouplingNumber
                              CALL List_ItemAdd(rankGlobalRowsLists(equationsIdx,rowRank)%ptr,rowListItem,err,error,*999)
                            ELSE
                              CALL FlagError("Global row is not owned by a domain.",err,error,*999)
                            ENDIF
                          ENDDO !globalRow
                        ELSE
                          CALL FlagError("Equations set boundary conditions is not associated.",err,error,*999)
                        ENDIF
                      ELSE
                        CALL FlagError("Equations set dependent field is not associated.",err,error,*999)
                      ENDIF
                    ELSE
                      CALL FlagError("Equations set row degree of freedom mappings is not associated.",err,error,*999)
                    ENDIF
                  ELSE
                    CALL FlagError("Equations mapping LHS mapping is not associated.",err,error,*999)
                  ENDIF
                ELSE
                  CALL FlagError("Equations equations mapping is not associated",err,error,*999)
                ENDIF
              ELSE
                CALL FlagError("Equations set equations is not associated.",err,error,*999)
              ENDIF
            ELSE
              CALL FlagError("Equations set is not associated.",err,error,*999)
            ENDIF
          ENDDO !equationsSetIdx
          !Now add in rows from any interface matrices
          DO interfaceConditionIdx=1,solverMapping%numberOfInterfaceConditions
            equationsIdx=equationsIdx+1
            interfaceCondition=>solverMapping%interfaceConditions(interfaceConditionIdx)%ptr
            IF(ASSOCIATED(interfaceCondition)) THEN
              SELECT CASE(interfaceCondition%method)
              CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD,INTERFACE_CONDITION_PENALTY_METHOD)
                interfaceEquations=>interfaceCondition%INTERFACE_EQUATIONS
                IF(ASSOCIATED(interfaceEquations)) THEN
                  interfaceMapping=>interfaceEquations%INTERFACE_MAPPING
                  IF(ASSOCIATED(interfaceMapping)) THEN
                    colDofsMapping=>interfaceMapping%COLUMN_DOFS_MAPPING
                    IF(ASSOCIATED(colDofsMapping)) THEN                    
                      lagrangeVariable=>interfaceMapping%LAGRANGE_VARIABLE
                      CALL BOUNDARY_CONDITIONS_VARIABLE_GET(boundaryConditions,lagrangeVariable,boundaryConditionsVariable, &
                        & err,error,*999)
                      IF(ASSOCIATED(boundaryConditionsVariable)) THEN
                        DO globalColumn=1,interfaceMapping%NUMBER_OF_GLOBAL_COLUMNS
                          !Find the rank that owns this global column
                          columnRank=-1
                          DO rankIdx=1,colDofsMapping%GLOBAL_TO_LOCAL_MAP(globalColumn)%NUMBER_OF_DOMAINS
                            IF(colDofsMapping%GLOBAL_TO_LOCAL_MAP(globalColumn)%LOCAL_TYPE(rankIdx)/=DOMAIN_LOCAL_GHOST) THEN
                              columnRank=colDofsMapping%GLOBAL_TO_LOCAL_MAP(globalColumn)%DOMAIN_NUMBER(rankIdx)
                              localColumn=colDofsMapping%GLOBAL_TO_LOCAL_MAP(globalColumn)%LOCAL_NUMBER(rankIdx)
                              EXIT
                            ENDIF
                          ENDDO !rankIdx
                          IF(columnRank>=0) THEN
                            globalDof=globalColumn
                            includeColumn=boundaryConditionsVariable%DOF_TYPES(globalDof)==BOUNDARY_CONDITION_DOF_FREE
                            rowListItem(1)=globalColumn
                            rowListItem(2)=localColumn
                            IF(includeColumn) THEN
                              rowListItem(3)=1
                              numberOfGlobalSolverRows=numberOfGlobalSolverRows+1
                              !Don't need to worry about ghosted rows.
                              IF(columnRank==myrank) numberOfLocalSolverRows=numberOfLocalSolverRows+1 !1-1 mapping
                            ELSE
                              rowListItem(3)=0                          
                            ENDIF !include column
                            rowListItem(4)=0
                            CALL List_ItemAdd(rankGlobalRowsLists(equationsIdx,columnRank)%ptr,rowListItem,err,error,*999)
                          ELSE
                            CALL FlagError("Global row is not owned by a domain.",err,error,*999)
                          ENDIF
                        ENDDO !globalColumn
                      ELSE
                        CALL FlagError("Boundary condition variable is not associated.",err,error,*999)
                      ENDIF
                    ELSE
                      CALL FlagError("Interface condition column degree of freedom mappings is not associated.",err,error,*999)
                    ENDIF
                  ELSE
                    CALL FlagError("Interface equations interface mapping is not associated.",err,error,*999)
                  ENDIF
                ELSE
                  CALL FlagError("Interface condition interface equations is not associated.",err,error,*999)
                ENDIF
              CASE(INTERFACE_CONDITION_AUGMENTED_LAGRANGE_METHOD)
                CALL FlagError("Not implemented.",err,error,*999)
              CASE(INTERFACE_CONDITION_POINT_TO_POINT_METHOD)
                CALL FlagError("Not implemented.",err,error,*999)
              CASE DEFAULT
                localError="The interface condition method of "//TRIM(NumberToVString(interfaceCondition%method,"*",err,error))// &
                  & " is invalid."
                CALL FlagError(localError,err,error,*999)
              END SELECT
            ELSE
              CALL FlagError("Interface condition is not associated.",err,error,*999)
            ENDIF
          ENDDO !interfaceConditionIdx

          !Sanity check.
          IF(numberOfLocalSolverRows==0) &
            & CALL FlagError("Invalid problem setup. The number of local solver rows is zero.",err,error,*999)
          IF(numberOfGlobalSolverRows==0) &
            & CALL FlagError("Invalid problem setup. The number of global solver rows is zero.",err,error,*999)

          !
          ! 3. We now know how many local and global rows are in the solver matrix. Loop over the rows in rank order and calculate
          !    the row mappings.
          !
          ! 3a Allocate and initialise the data structures
          !

          !Allocate memory for the rows mapping
          !Allocate the solver rows to equations set maps
          ALLOCATE(solverMapping%solverRowToEquationsRowsMap(numberOfLocalSolverRows),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate solver mapping solver row to equation rows map.",err,error,*999)
          !Set the number of rows
          solverMapping%numberOfRows=numberOfLocalSolverRows
          solverMapping%numberOfGlobalRows=numberOfGlobalSolverRows
          !Allocate the solver rows domain mapping
          ALLOCATE(solverMapping%rowDofsMapping,STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate solver mapping row dofs mapping.",err,error,*999)
!!TODO: what is the real number of domains for a solver???
          CALL DOMAIN_MAPPINGS_MAPPING_INITIALISE(solverMapping%rowDofsMapping,COMPUTATIONAL_ENVIRONMENT% &
            & NUMBER_COMPUTATIONAL_NODES,err,error,*999)
          rowDomainMapping=>solverMapping%rowDofsMapping
          ALLOCATE(rowDomainMapping%GLOBAL_TO_LOCAL_MAP(numberOfGlobalSolverRows),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate row dofs mapping global to local map.",err,error,*999)
          rowDomainMapping%NUMBER_OF_GLOBAL=numberOfGlobalSolverRows

          !Initialise the equations sets to solver maps
          DO equationsSetIdx=1,solverMapping%numberOfEquationsSets

            !Note that pointers have been checked for association above
            equationsSet=>solverMapping%equationsSets(equationsSetIdx)%ptr
            equations=>equationsSet%equations
            equationsMapping=>equations%EQUATIONS_MAPPING
            
            !Allocate the equations set to solver maps for solver matrix (sm) indexing
            ALLOCATE(solverMapping%equationsSetToSolverMap(equationsSetIdx)%equationsToSolverMatrixMapsSm( &
              & solverMapping%numberOfSolverMatrices),STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate equations set to solver map equations to solver matrix maps sm.", &
              & err,error,*999)
            DO solverMatrixIdx=1,solverMapping%numberOfSolverMatrices
              CALL SolverMapping_EquationsToSolverMatrixMapsSmInitialise(solverMapping%equationsSetToSolverMap( &
                & equationsSetIdx)%equationsToSolverMatrixMapsSm(solverMatrixIdx),err,error,*999)
            ENDDO !solverMatrixIdx
            
            !Allocate the equations row to solver rows maps
            ALLOCATE(solverMapping%equationsSetToSolverMap(equationsSetIdx)%equationsRowToSolverRowsMaps( &
              & equationsMapping%TOTAL_NUMBER_OF_ROWS),STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate equations set to solver map equations row to solver rows maps.", &
              & err,error,*999)
            DO equationsRowNumber=1,equationsMapping%TOTAL_NUMBER_OF_ROWS
              !Initialise
              CALL SolverMapping_EquationsRowToSolverRowsMapInitialise(solverMapping%equationsSetToSolverMap( &
                & equationsSetIdx)%equationsRowToSolverRowsMaps(equationsRowNumber),err,error,*999)
            ENDDO
            
          ENDDO !equationsSetIdx
          
          !Initialise the interface condition to solver maps
          DO interfaceConditionIdx=1,solverMapping%numberOfInterfaceConditions
            
            !Note that pointers have been checked for association above
            interfaceCondition=>solverMapping%interfaceConditions(interfaceConditionIdx)%ptr
            interfaceEquations=>interfaceCondition%INTERFACE_EQUATIONS
            interfaceMapping=>interfaceEquations%INTERFACE_MAPPING
            
            !Allocate the interface to solver maps for solver matrix (sm) indexing
            ALLOCATE(solverMapping%interfaceConditionToSolverMap(interfaceConditionIdx)%interfaceToSolverMatrixMapsSm( &
              & solverMapping%numberOfSolverMatrices),STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate interface to solver map interface to solver matrix maps sm.", &
              & err,error,*999)
            DO solverMatrixIdx=1,solverMapping%numberOfSolverMatrices
              CALL SolverMapping_InterfaceToSolverMatrixMapsSmInitialise(solverMapping%interfaceConditionToSolverMap( &
                & interfaceConditionIdx)%interfaceToSolverMatrixMapsSm(solverMatrixIdx),err,error,*999)
            ENDDO !solverMatrixIdx
            
            !Allocate the interface to solver maps for interface matrix (im) indexing
            ALLOCATE(solverMapping%interfaceConditionToSolverMap(interfaceConditionIdx)% &
              & interfaceToSolverMatrixMapsIm(interfaceMapping%NUMBER_OF_INTERFACE_MATRICES),STAT=err)
            IF(err/=0) &
              & CALL FlagError("Could not allocate interface to solver map equations to solver matrix maps im.", &
              & err,error,*999)
            DO interfaceMatrixIdx=1,interfaceMapping%NUMBER_OF_INTERFACE_MATRICES
              CALL SolverMapping_InterfaceToSolverMatrixMapsImInitialise(solverMapping%interfaceConditionToSolverMap( &
                & interfaceConditionIdx)%interfaceToSolverMatrixMapsIm(interfaceMatrixIdx),err,error,*999)
                      
              !Allocate the interfafce row to solver row maps
              ALLOCATE(solverMapping%interfaceConditionToSolverMap(interfaceConditionIdx)% &
                & interfaceToSolverMatrixMapsIm(interfaceMatrixIdx)%interfaceRowToSolverRowsMap( &
                & interfaceMapping%INTERFACE_MATRIX_ROWS_TO_VAR_MAPS(interfaceMatrixIdx)%TOTAL_NUMBER_OF_ROWS),STAT=err)
              IF(err/=0) CALL FlagError("Could not allocate interface condition to solver map interface row to solver row map.", &
                & err,error,*999)
              DO interfaceRowNumber=1,interfaceMapping%INTERFACE_MATRIX_ROWS_TO_VAR_MAPS(interfaceMatrixIdx)%TOTAL_NUMBER_OF_ROWS
                CALL SolverMapping_InterfaceRowToSolverRowsMapInitialise(solverMapping%interfaceConditionToSolverMap( &
                  & interfaceConditionIdx)%interfaceToSolverMatrixMapsIm(interfaceMatrixIdx)% &
                  & interfaceRowToSolverRowsMap(interfaceRowNumber),err,error,*999)                
              ENDDO !interfaceRowNumber
              
            ENDDO !interfaceMatrixIdx

            !Allocate the interface column to solver row maps
            ALLOCATE(solverMapping%interfaceConditionToSolverMap(interfaceConditionIdx)% &
              & interfaceColumnToSolverRowsMaps(interfaceMapping%TOTAL_NUMBER_OF_COLUMNS),STAT=err)
            IF(err/=0)  &
              & CALL FlagError("Could not allocate interface condition to solver map interface column to solver row map.", &
              & err,error,*999)
            DO interfaceColNumber=1,interfaceMapping%TOTAL_NUMBER_OF_COLUMNS
              !Initialise
              CALL SolverMapping_InterfaceColToSolverRowsMapInitialise(solverMapping%interfaceConditionToSolverMap( &
                & interfaceConditionIdx)%interfaceColumnToSolverRowsMaps(interfaceColNumber),err,error,*999)
            ENDDO !interfaceColNumber
            
          ENDDO !interfaceConditionIdx

          !
          ! 3b Now calculate the mappings for each solver row <-> equations row & interface row/column in rank order.
          !
          
          !Calculate the row mappings
          numberOfGlobalSolverRows=0
          !Make a "dof coupling" for rows that aren't coupled
          ALLOCATE(dummyDofCoupling%globalDofs(1),stat=err)
          IF(err/=0) CALL FlagError("Could not allocate dummy DOF coupling DOFs.",err,error,*999)
          ALLOCATE(dummyDofCoupling%localDofs(1),stat=err)
          IF(err/=0) CALL FlagError("Could not allocate dummy DOF coupling local Dofs.",err,error,*999)
          ALLOCATE(dummyDofCoupling%coefficients(1),stat=err)
          IF(err/=0) CALL FlagError("Could not allocate dummy DOF coupling values.",err,error,*999)
          dummyDofCoupling%numberOfDofs=1
          !Loop over the ranks to  ensure that the lowest ranks have the lowest numbered solver variables
          DO rank=0,COMPUTATIONAL_ENVIRONMENT%NUMBER_COMPUTATIONAL_NODES-1
            numberOfLocalSolverRows=0
            equationsIdx=0
            DO equationsSetIdx=1,solverMapping%numberOfEquationsSets
              equationsIdx=equationsIdx+1

              !Get rows list
              CALL List_Sort(rankGlobalRowsLists(equationsIdx,rank)%ptr,err,error,*999)
              CALL List_DetachAndDestroy(rankGlobalRowsLists(equationsIdx,rank)%ptr,numberOfRankRows,rankGlobalRowsList, &
                & err,error,*999)

              !Note that pointers have been checked for association above
              equationsSet=>solverMapping%equationsSets(equationsSetIdx)%ptr
              equations=>equationsSet%equations
              equationsMapping=>equations%EQUATIONS_MAPPING
              dynamicMapping=>equationsMapping%DYNAMIC_MAPPING
              linearMapping=>equationsMapping%LINEAR_MAPPING
              nonlinearMapping=>equationsMapping%NONLINEAR_MAPPING

              !Loop over the global rows for this rank.
              DO globalRowIdx=1,numberOfRankRows
                globalRow=rankGlobalRowsList(1,globalRowIdx)
                localRow=rankGlobalRowsList(2,globalRowIdx)
                includeRow=rankGlobalRowsList(3,globalRowIdx)==1
                constrainedDof=rankGlobalRowsList(3,globalRowIdx)==2
                globalDofCouplingNumber=rankGlobalRowsList(4,globalRowIdx)
                IF(globalDofCouplingNumber>0) THEN
                  rowEquationRows=>rowCouplings%dofCouplings(globalDofCouplingNumber)%ptr
                  IF(ASSOCIATED(rowEquationRows)) THEN
                    numberRowEquationsRows=rowEquationRows%numberOfDofs
                  ELSE
                    CALL FlagError("Dof coupling is not associated for global dof coupling number "// &
                      & TRIM(NumberToVstring(globalDofCouplingNumber,"*",err,error))//".",err,error,*999)
                  END IF
                ELSE
                  numberRowEquationsRows=1
                  dummyDofCoupling%globalDofs(1)=globalRow
                  dummyDofCoupling%localDofs(1)=localRow
                  dummyDofCoupling%coefficients(1)=1.0_DP
                  rowEquationRows=>dummyDofCoupling
                END IF
                IF(includeRow) THEN
                  numberOfGlobalSolverRows=numberOfGlobalSolverRows+1
                  numberOfLocalSolverRows=numberOfLocalSolverRows+1
                  !Set up the row domain mappings.
                  !There are no ghosted rows for the solver matrices so there is only one domain for the global to local map.
                  !Initialise
                  CALL DOMAIN_MAPPINGS_MAPPING_GLOBAL_INITIALISE(rowDomainMapping%GLOBAL_TO_LOCAL_MAP(numberOfGlobalSolverRows), &
                    & err,error,*999)
                  !Allocate the global to local map arrays
                  ALLOCATE(rowDomainMapping%GLOBAL_TO_LOCAL_MAP(numberOfGlobalSolverRows)%LOCAL_NUMBER(1),STAT=err)
                  IF(err/=0) CALL FlagError("Could not allocate row global to local map local number.",err,error,*999)
                  ALLOCATE(rowDomainMapping%GLOBAL_TO_LOCAL_MAP(numberOfGlobalSolverRows)%DOMAIN_NUMBER(1),STAT=err)
                  IF(err/=0) CALL FlagError("Could not allocate row global to local map domain number.",err,error,*999)
                  ALLOCATE(rowDomainMapping%GLOBAL_TO_LOCAL_MAP(numberOfGlobalSolverRows)%LOCAL_TYPE(1),STAT=err)
                  IF(err/=0) CALL FlagError("Could not allocate row global to local map local type.",err,error,*999)
                  !Set the global to local mappings
                  rowDomainMapping%GLOBAL_TO_LOCAL_MAP(numberOfGlobalSolverRows)%NUMBER_OF_DOMAINS=1
                  rowDomainMapping%GLOBAL_TO_LOCAL_MAP(numberOfGlobalSolverRows)%LOCAL_NUMBER(1)=numberOfLocalSolverRows
                  rowDomainMapping%GLOBAL_TO_LOCAL_MAP(numberOfGlobalSolverRows)%DOMAIN_NUMBER(1)=rank
                  rowDomainMapping%GLOBAL_TO_LOCAL_MAP(numberOfGlobalSolverRows)%LOCAL_TYPE(1)=DOMAIN_LOCAL_INTERNAL
                  IF(rank==myrank) THEN
                    !If this is my rank then set up the solver->equations and equations->solver row mappings
                    
                    !Set up the solver row -> equations row mappings. Will need to look
                    !At the interface conditions for this equations set later.
                    !Initialise
                    CALL SolverMapping_SolverRowToEquationsMapsInitialise(solverMapping%solverRowToEquationsRowsMap( &
                      & numberOfLocalSolverRows),err,error,*999)

                    ALLOCATE(solverMapping%solverRowToEquationsRowsMap(numberOfLocalSolverRows)% &
                      & equationsIndex(numberRowEquationsRows),STAT=err)
                    IF(err/=0) CALL FlagError("Could not allocate solver row to equations rows equations index.",err,error,*999)
                    ALLOCATE(solverMapping%solverRowToEquationsRowsMap(numberOfLocalSolverRows)% &
                      & rowcolNumber(numberRowEquationsRows),STAT=err)
                    IF(err/=0) &
                      & CALL FlagError("Could not allocate solver row to equations rows row/col number.",err,error,*999)
                    ALLOCATE(solverMapping%solverRowToEquationsRowsMap(numberOfLocalSolverRows)% &
                      & couplingCoefficients(numberRowEquationsRows),STAT=err)
                    IF(err/=0) CALL FlagError("Could not allocate solver row to equations rows coupling coefficients.", &
                      & err,error,*999)
                    !Set the mappings for the first equations DOF, the rest will be set up later using the DOF constraints
                    solverMapping%solverRowToEquationsRowsMap(numberOfLocalSolverRows)%numberOfEquationsSets=numberRowEquationsRows
                    DO rowEquationsRowIdx=1,numberRowEquationsRows
                      solverMapping%solverRowToEquationsRowsMap(numberOfLocalSolverRows)%equationsIndex(rowEquationsRowIdx)= &
                        & equationsSetIdx
                      solverMapping%solverRowToEquationsRowsMap(numberOfLocalSolverRows)%rowcolNumber(rowEquationsRowIdx)= &
                        rowEquationRows%localDofs(rowEquationsRowIdx)
                      solverMapping%solverRowToEquationsRowsMap(numberOfLocalSolverRows)%couplingCoefficients(rowEquationsRowIdx)= &
                        & rowEquationRows%coefficients(rowEquationsRowIdx)
                    ENDDO !rowEquationsRowIdx
                    !Set up the equations row -> solver row mappings
                    DO rowEquationsRowIdx=1,numberRowEquationsRows
                      equationsRow=rowEquationRows%localDofs(rowEquationsRowIdx)
                      !Allocate the equations row to solver row mappings arrays
                      IF(ALLOCATED(solverMapping%equationsSetToSolverMap(equationsSetIdx)% &
                          & equationsRowToSolverRowsMaps(equationsRow)%solverRows)) THEN
                        CALL FlagError("Equations row to solver row map is already allocated, "// &
                          & "mapping an equations row to multiple solver rows is not yet implemented.",err,error,*999)
                      END IF
                      ALLOCATE(solverMapping%equationsSetToSolverMap(equationsSetIdx)%equationsRowToSolverRowsMaps( &
                        & localRow)%solverRows(1),STAT=err)
                      IF(err/=0) CALL FlagError("Could not allocate equations row to solver rows maps solver rows.", &
                        & err,error,*999)
                      ALLOCATE(solverMapping%equationsSetToSolverMap(equationsSetIdx)%equationsRowToSolverRowsMaps( &
                        & localRow)%couplingCoefficients(1),STAT=err)
                      IF(err/=0) CALL FlagError("Could not allocate equations row to solver rows maps solver rows.", &
                        & err,error,*999)
                      !Set the mappings
                      solverMapping%equationsSetToSolverMap(equationsSetIdx)%equationsRowToSolverRowsMaps( &
                        & equationsRow)%numberOfSolverRows=1
                      solverMapping%equationsSetToSolverMap(equationsSetIdx)%equationsRowToSolverRowsMaps( &
                        & equationsRow)%solverRows(1)=numberOfLocalSolverRows
                      solverMapping%equationsSetToSolverMap(equationsSetIdx)%equationsRowToSolverRowsMaps( &
                        & equationsRow)%couplingCoefficients(1)=rowEquationRows%coefficients(rowEquationsRowIdx)
                      !Now set up any interface condition rows to solver rows that affect this equations set.
                      DO interfaceConditionIdx=1,solverMapping%equationsSetToSolverMap(equationsSetIdx)%numberOfInterfaceConditions
                        interfaceConditionIdx2=solverMapping%equationsSetToSolverMap(equationsSetIdx)% &
                          & equationsToSolverMatrixMapsInterface(interfaceConditionIdx)%interfaceConditionIndex
                        interfaceMatrixIdx=solverMapping%equationsSetToSolverMap(equationsSetIdx)% &
                          & equationsToSolverMatrixMapsInterface(interfaceConditionIdx)%interfaceMatrixNumber
                        !Set the mappings
                        solverMapping%interfaceConditionToSolverMap(interfaceConditionIdx2)% &
                          & interfaceToSolverMatrixMapsIm(interfaceMatrixIdx)%interfaceRowToSolverRowsMap( &
                          & equationsRow)%numberOfSolverRows=1
                        solverMapping%interfaceConditionToSolverMap(interfaceConditionIdx2)% &
                          & interfaceToSolverMatrixMapsIm(interfaceMatrixIdx)%interfaceRowToSolverRowsMap( &
                          & EquationsRow)%solverRow=numberOfLocalSolverRows
                        solverMapping%interfaceConditionToSolverMap(interfaceConditionIdx2)% &
                          & interfaceToSolverMatrixMapsIm(interfaceMatrixIdx)%interfaceRowToSolverRowsMap( &
                          & EquationsRow)%couplingCoefficient=rowEquationRows%coefficients(rowEquationsRowIdx)              
                      ENDDO !interfaceConditionIdx
                    ENDDO !rowEquationsRowIdx
                  ENDIF !rank==my rank
                ELSE IF(constrainedDof) THEN
                  !Do nothing, the row mapping for this equations row will be set up with
                  !the first equations row mapped to the solver row
                ELSE 
                  !Set up the equations row -> solver row mappings. 
                  !Note that in this case these mappings are set to zero since these equation rows don't appear in the solver matrices
                  IF(rank==myrank) THEN
                    !Set the mappings
                    solverMapping%equationsSetToSolverMap(equationsSetIdx)%equationsRowToSolverRowsMaps(localRow)% &
                      & numberOfSolverRows=0
                    !Now set up any interface condition rows to solver rows that affect this equations set.
                    DO interfaceConditionIdx=1,solverMapping%equationsSetToSolverMap(equationsSetIdx)%numberOfInterfaceConditions
                      interfaceConditionIdx2=solverMapping%equationsSetToSolverMap(equationsSetIdx)% &
                        & equationsToSolverMatrixMapsInterface(interfaceConditionIdx)%interfaceConditionIndex
                      interfaceMatrixIdx=solverMapping%equationsSetToSolverMap(equationsSetIdx)% &
                        & equationsToSolverMatrixMapsInterface(interfaceConditionIdx)%interfaceMatrixNumber
                      !Set the mappings
                      solverMapping%interfaceConditionToSolverMap(interfaceConditionIdx2)% &
                        & interfaceToSolverMatrixMapsIm(interfaceMatrixIdx)%interfaceRowToSolverRowsMap( &
                        & localRow)%numberOfSolverRows=0
                    ENDDO !interfaceConditionIdx
                  ENDIF !rank==my rank
                ENDIF !include row
              ENDDO !globalRowIdx
              IF(ALLOCATED(rankGlobalRowsList)) DEALLOCATE(rankGlobalRowsList)
            ENDDO !equationsSetIdx
            DO interfaceConditionIdx=1,solverMapping%numberOfInterfaceConditions
              equationsIdx=equationsIdx+1

              !Get rows list
              CALL List_Sort(rankGlobalRowsLists(equationsIdx,rank)%ptr,err,error,*999)
              CALL List_DetachAndDestroy(rankGlobalRowsLists(equationsIdx,rank)%ptr,numberOfRankRows,rankGlobalRowsList, &
                & err,error,*999)

              !Note that pointers have been checked for association above
              interfaceCondition=>solverMapping%interfaceConditions(interfaceConditionIdx)%ptr
              interfaceEquations=>interfaceCondition%INTERFACE_EQUATIONS
              interfaceMapping=>interfaceEquations%INTERFACE_MAPPING

              !Loop over the global rows for this rank.
              DO globalRowIdx=1,numberOfRankRows
                globalColumn=rankGlobalRowsList(1,globalRowIdx)
                localColumn=rankGlobalRowsList(2,globalRowIdx)
                includeColumn=rankGlobalRowsList(3,globalRowIdx)==1
                constrainedDof=rankGlobalRowsList(3,globalRowIdx)==2
                IF(includeColumn) THEN
                  numberOfGlobalSolverRows=numberOfGlobalSolverRows+1
                  numberOfLocalSolverRows=numberOfLocalSolverRows+1
                  !Set up the row domain mappings.
                  !There are no ghosted rows for the solver matrices so there is only one domain for the global to local map.
                  !Initialise
                  CALL DOMAIN_MAPPINGS_MAPPING_GLOBAL_INITIALISE(rowDomainMapping%GLOBAL_TO_LOCAL_MAP(numberOfGlobalSolverRows), &
                    & err,error,*999)
                  !Allocate the global to local map arrays
                  ALLOCATE(rowDomainMapping%GLOBAL_TO_LOCAL_MAP(numberOfGlobalSolverRows)%LOCAL_NUMBER(1),STAT=err)
                  IF(err/=0) CALL FlagError("Could not allocate row global to local map local number.",err,error,*999)
                  ALLOCATE(rowDomainMapping%GLOBAL_TO_LOCAL_MAP(numberOfGlobalSolverRows)%DOMAIN_NUMBER(1),STAT=err)
                  IF(err/=0) CALL FlagError("Could not allocate row global to local map domain number.",err,error,*999)
                  ALLOCATE(rowDomainMapping%GLOBAL_TO_LOCAL_MAP(numberOfGlobalSolverRows)%LOCAL_TYPE(1),STAT=err)
                  IF(err/=0) CALL FlagError("Could not allocate row global to local map local type.",err,error,*999)
                  !Set the global to local mappings
                  rowDomainMapping%GLOBAL_TO_LOCAL_MAP(numberOfGlobalSolverRows)%NUMBER_OF_DOMAINS=1
                  rowDomainMapping%GLOBAL_TO_LOCAL_MAP(numberOfGlobalSolverRows)%LOCAL_NUMBER(1)= &
                    & numberOfLocalSolverRows
                  rowDomainMapping%GLOBAL_TO_LOCAL_MAP(numberOfGlobalSolverRows)%DOMAIN_NUMBER(1)=rank
                  rowDomainMapping%GLOBAL_TO_LOCAL_MAP(numberOfGlobalSolverRows)%LOCAL_TYPE(1)=DOMAIN_LOCAL_INTERNAL
                  !If this is my rank then set up the solver->equations and equations->solver row mappings
                  IF(rank==myrank) THEN
                    !Set the interface column/row -> solver row mappings
                    !Note that for populating solverMapping%solverRowToEquationsRowsMap(numberOfLocalSolverRows)%rowColNumber(i)
                    !If the row are equations set rows this is the i'th row number that the solver row is mapped to.
                    !If the rows are interface rows (which is the case here) then this is the i'th column number that the solver
                    !row is mapped to.
                    !Initialise
                    CALL SolverMapping_SolverRowToEquationsMapsInitialise(solverMapping%solverRowToEquationsRowsMap( &
                      & numberOfLocalSolverRows),err,error,*999)

                    ALLOCATE(solverMapping%solverRowToEquationsRowsMap(numberOfLocalSolverRows)%rowcolNumber(1),STAT=err)
                    IF(err/=0) &
                      & CALL FlagError("Could not allocate solver row to equations rows row/col number.",err,error,*999)
                    ALLOCATE(solverMapping%solverRowToEquationsRowsMap(numberOfLocalSolverRows)%couplingCoefficients(1), &
                      & STAT=err)
                    IF(err/=0) CALL FlagError("Could not allocate solver row to equations rows coupling coefficients.", &
                      & err,error,*999)
                    !Set the mappings
                    solverMapping%solverRowToEquationsRowsMap(numberOfLocalSolverRows)%interfaceConditionIndex= &
                      & interfaceConditionIdx
                    solverMapping%solverRowToEquationsRowsMap(numberOfLocalSolverRows)%rowcolNumber(1)=localColumn
                    solverMapping%solverRowToEquationsRowsMap(numberOfLocalSolverRows)%couplingCoefficients(1)=1.0_DP
                    !Set up the interface col -> solver row mappings
                    solverMapping%interfaceConditionToSolverMap(interfaceConditionIdx)% &
                      & interfaceColumnToSolverRowsMaps(localColumn)%numberOfSolverRows=1
                    solverMapping%interfaceConditionToSolverMap(interfaceConditionIdx)% &
                      & interfaceColumnToSolverRowsMaps(localColumn)%solverRow=numberOfLocalSolverRows
                    solverMapping%interfaceConditionToSolverMap(interfaceConditionIdx)% &
                      & interfaceColumnToSolverRowsMaps(localColumn)%couplingCoefficient=1.0_DP
                  ENDIF !rank==my rank
                ELSE IF(constrainedDof) THEN
                  CALL FlagError("Constrained DOFs have not been implemented for Lagrange variables.",err,error,*999)
                ELSE 
                  !Set the interface column/row -> solver row mappings
                  IF(rank==myrank) THEN
                    !Set the interface column -> solver row mappings
                    solverMapping%interfaceConditionToSolverMap(interfaceConditionIdx)% &
                      & interfaceColumnToSolverRowsMaps(localColumn)%numberOfSolverRows=0
                    SELECT CASE(interfaceCondition%method)
                    CASE(INTERFACE_CONDITION_PENALTY_METHOD)
                      !Set up the solver row <-> interface row mappings
                      !Penalty matrix is the last interface matrix
                      interfaceMatrixIdx=interfaceMapping%NUMBER_OF_INTERFACE_MATRICES
                      !Set the mappings
                      solverMapping%interfaceConditionToSolverMap(interfaceConditionIdx)% &
                        & interfaceToSolverMatrixMapsIm(interfaceMatrixIdx)%interfaceRowToSolverRowsMap(localColumn)% &
                        & numberOfSolverRows=0
                    ENDSELECT
                  ENDIF
                ENDIF
              ENDDO !globalRowIdx
              IF(ALLOCATED(rankGlobalRowsList)) DEALLOCATE(rankGlobalRowsList)              
            ENDDO !interfaceConditionIdx            
          ENDDO !rank
          CALL SolverDofCouplings_Finalise(rowCouplings,err,error,*999)
          
          CALL DOMAIN_MAPPINGS_LOCAL_FROM_GLOBAL_CALCULATE(rowDomainMapping,err,error,*999)
          IF(ALLOCATED(variablesList)) DEALLOCATE(variablesList)
          IF(ALLOCATED(variableProcessed)) DEALLOCATE(variableProcessed)
          IF(ALLOCATED(numberOfVariableGlobalSolverDofs)) DEALLOCATE(numberOfVariableGlobalSolverDofs)
          IF(ALLOCATED(numberOfVariableLocalSolverDofs)) DEALLOCATE(numberOfVariableLocalSolverDofs)
          IF(ALLOCATED(totalNumberOfVariableLocalSolverDofs)) DEALLOCATE(totalNumberOfVariableLocalSolverDofs)
          
          !
          !--- Column mappings ---
          !
          ! 4. Calculate the number of local and global columns in the solver matrix. Do this by calculating the list of columns
          !    for each rank so that we can determine the column numbers in rank order.
          !
          !Allocate solver column to equations sets mapping array
          ALLOCATE(solverMapping%solverColToEquationsColsMap(solverMapping%numberOfSolverMatrices),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate solver mapping solver column to equations column maps.",err,error,*999)
          !Allocate the column variable lists
          ALLOCATE(solverMapping%columnVariablesList(solverMapping%numberOfSolverMatrices),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate solver mapping column variables list.",err,error,*999)
          
          CALL SolverDofCouplings_Initialise(columnCouplings,err,error,*999)
         
          !Calculate the column mappings for each solver matrix
          DO solverMatrixIdx=1,solverMapping%numberOfSolverMatrices

            !Initialise the column variables list
            CALL SolverMapping_VariablesInitialise(solverMapping%columnVariablesList(solverMatrixIdx),err,error,*999)
            !
            ! 4a Calculate the list of field variables involved in the columns of the solver matrix
            !
            !Compute the order of variables for the solver matrices          
            CALL List_DetachAndDestroy(solverMapping%createValuesCache%equationsColVariablesList(solverMatrixIdx)%ptr, &
              & numberOfEquationsVariables,equationsVariables,err,error,*999)
            CALL List_DetachAndDestroy(solverMapping%createValuesCache%interfaceColVariablesList(solverMatrixIdx)%ptr, &
              & numberOfInterfaceVariables,interfaceVariables,err,error,*999)
            ALLOCATE(solverMapping%columnVariablesList(solverMatrixIdx)%variables(numberOfEquationsVariables+ &
              & numberOfInterfaceVariables),STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate column variables list variables.",err,error,*999)
            solverMapping%columnVariablesList(solverMatrixIdx)%numberOfVariables=numberOfEquationsVariables+ &
              & numberOfInterfaceVariables
            ALLOCATE(variablesList(numberOfEquationsVariables+numberOfInterfaceVariables),STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate variables list.",err,error,*999)
            ALLOCATE(variableProcessed(numberOfEquationsVariables+numberOfInterfaceVariables),STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate variable processed.",err,error,*999)
            variableProcessed=.FALSE.
            solverVariableIdx=0
            DO variableIdx=1,numberOfEquationsVariables
              solverVariableIdx=solverVariableIdx+1
              CALL SolverMapping_VariableInitialise(solverMapping%columnVariablesList(solverMatrixIdx)% &
                & variables(solverVariableIdx),err,error,*999)
              equationsSetIdx=equationsVariables(1,variableIdx)
              equationsSet=>solverMapping%equationsSets(equationsSetIdx)%ptr
              IF(ASSOCIATED(equationsSet)) THEN
                dependentField=>equationsSet%dependent%DEPENDENT_FIELD
                IF(ASSOCIATED(dependentField)) THEN
                  variableType=equationsVariables(2,variableIdx)
                  variable=>dependentField%VARIABLE_TYPE_MAP(variableType)%ptr
                  IF(ASSOCIATED(variable)) THEN
                    solverMapping%columnVariablesList(solverMatrixIdx)%variables(solverVariableIdx)%variable=>variable
                    solverMapping%columnVariablesList(solverMatrixIdx)%variables(solverVariableIdx)%variableType=variableType
                    NULLIFY(variablesList(solverVariableIdx)%ptr)
                    CALL List_CreateStart(variablesList(solverVariableIdx)%ptr,err,error,*999)
                    CALL List_DataTypeSet(variablesList(solverVariableIdx)%ptr,LIST_INTG_TYPE,err,error,*999)
                    CALL List_DataDimensionSet(variablesList(solverVariableIdx)%ptr,3,err,error,*999)
                    CALL List_KeyDimensionSet(variablesList(solverVariableIdx)%ptr,1,err,error,*999)
                    CALL List_CreateFinish(variablesList(solverVariableIdx)%ptr,err,error,*999)
                  ELSE
                    CALL FlagError("Dependent field variable is not associated.",err,error,*999)
                  ENDIF
                ELSE
                  CALL FlagError("Equations set dependent field is not associated.",err,error,*999)
                ENDIF
              ELSE
                CALL FlagError("Equations set is not associated.",err,error,*999)
              ENDIF
            ENDDO !variableIdx
            IF(ALLOCATED(equationsVariables)) DEALLOCATE(equationsVariables)
            DO variableIdx=1,numberOfInterfaceVariables
              solverVariableIdx=solverVariableIdx+1
              CALL SolverMapping_VariableInitialise(solverMapping%columnVariablesList(solverMatrixIdx)% &
                & variables(solverVariableIdx),err,error,*999)
              interfaceConditionIdx=interfaceVariables(1,variableIdx)
              interfaceCondition=>solverMapping%interfaceConditions(interfaceConditionIdx)%ptr
              IF(ASSOCIATED(interfaceCondition)) THEN
                IF(ASSOCIATED(interfaceCondition%lagrange)) THEN
                  lagrangeField=>interfaceCondition%lagrange%LAGRANGE_FIELD
                  IF(ASSOCIATED(lagrangeField)) THEN
                    variableType=interfaceVariables(2,variableIdx)
                    variable=>lagrangeField%VARIABLE_TYPE_MAP(variableType)%ptr
                    IF(ASSOCIATED(variable)) THEN
                      solverMapping%columnVariablesList(solverMatrixIdx)%variables(solverVariableIdx)%variable=>variable
                      solverMapping%columnVariablesList(solverMatrixIdx)%variables(solverVariableIdx)%variableType=variableType
                      NULLIFY(variablesList(solverVariableIdx)%ptr)
                      CALL List_CreateStart(variablesList(solverVariableIdx)%ptr,err,error,*999)
                      CALL List_DataTypeSet(variablesList(solverVariableIdx)%ptr,LIST_INTG_TYPE,err,error,*999)
                      CALL List_DataDimensionSet(variablesList(solverVariableIdx)%ptr,3,err,error,*999)
                      CALL List_KeyDimensionSet(variablesList(solverVariableIdx)%ptr,1,err,error,*999)
                      CALL List_CreateFinish(variablesList(solverVariableIdx)%ptr,err,error,*999)
                    ELSE
                      CALL FlagError("Lagrange field variable is not associated.",err,error,*999)
                    ENDIF
                  ELSE
                    CALL FlagError("Interface condition Lagrange field is not associated.",err,error,*999)
                  ENDIF
                ELSE
                  CALL FlagError("Interface condition Lagrange is not associated.",err,error,*999)
                ENDIF
              ELSE
                CALL FlagError("Interface condition is not associated.",err,error,*999)
              ENDIF
            ENDDO !variableIdx
            IF(ALLOCATED(interfaceVariables)) DEALLOCATE(interfaceVariables)
            
            !Initialise solver column to equations sets mapping array
            CALL SolverMapping_SolverColToEquationsMapsInitialise(solverMapping%solverColToEquationsColsMap( &
              & solverMatrixIdx),err,error,*999)
            solverMapping%solverColToEquationsColsMap(solverMatrixIdx)%solverMatrixNumber=solverMatrixIdx
            solverMapping%solverColToEquationsColsMap(solverMatrixIdx)%solverMapping=>solverMapping
            !Allocate the solver col to equations set maps array
            ALLOCATE(solverMapping%solverColToEquationsColsMap(solverMatrixIdx)%solverColToEquationsSetMaps( &
              & solverMapping%numberOfEquationsSets),STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate solver col to equations map solver col to equation set maps.", &
              & err,error,*999)
            !Allocate the solver col to interface maps array
            ALLOCATE(solverMapping%solverColToEquationsColsMap(solverMatrixIdx)%solverColToInterfaceMaps( &
              & solverMapping%numberOfInterfaceConditions),STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate solver col to equations map solver col to interface maps.", &
              & err,error,*999)
            !Presort the column numbers by rank.
            !RANK_GLOBAL_COLS_LISTS(dof_type, equations_idx, variable_idx, rank_idx)
            !dof_type is 1 for domain local DOFs and 2 for ghost DOFs
            ALLOCATE(rankGlobalColsLists(2,solverMapping%numberOfEquationsSets+solverMapping%numberOfInterfaceConditions, &
              & solverMapping%columnVariablesList(solverMatrixIdx)%numberOfVariables, &
              & 0:COMPUTATIONAL_ENVIRONMENT%NUMBER_COMPUTATIONAL_NODES-1),STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate rank global columns lists.",err,error,*999)
            DO rank=0,COMPUTATIONAL_ENVIRONMENT%NUMBER_COMPUTATIONAL_NODES-1
              DO solverVariableIdx=1,solverMapping%columnVariablesList(solverMatrixIdx)%numberOfVariables
                DO equationsIdx=1,solverMapping%numberOfEquationsSets+solverMapping%numberOfInterfaceConditions
                  DO dofType=1,2
                    NULLIFY(rankGlobalColsLists(dofType,equationsIdx,solverVariableIdx,rank)%ptr)
                    CALL List_CreateStart(rankGlobalColsLists(dofType,equationsIdx,solverVariableIdx,rank)%ptr,err,error,*999)
                    CALL List_DataTypeSet(rankGlobalColsLists(dofType,equationsIdx,solverVariableIdx,rank)%ptr, &
                      & LIST_INTG_TYPE,err,error,*999)
                    !Set approximate size for the number of columns per variable.
                    CALL List_InitialSizeSet(rankGlobalColsLists(dofType,equationsIdx,solverVariableIdx,rank)%ptr, &
                      & solverMapping%columnVariablesList(solverMatrixIdx)%variables(solverVariableIdx)%variable% &
                      & TOTAL_NUMBER_OF_DOFS,err,error,*999)
                    CALL List_DataDimensionSet(rankGlobalColsLists(dofType,equationsIdx,solverVariableIdx,rank)%ptr,5, &
                      & err,error,*999)
                    CALL List_KeyDimensionSet(rankGlobalColsLists(dofType,equationsIdx,solverVariableIdx,rank)%ptr,1, &
                      & err,error,*999)
                    CALL List_CreateFinish(rankGlobalColsLists(dofType,equationsIdx,solverVariableIdx,rank)%ptr, &
                      & err,error,*999)
                  ENDDO !dofType
                ENDDO !equationsIdx
              ENDDO !solverVariableIdx
            ENDDO !rank

            !Allocate sub-matrix information
            !SUB_MATRIX_INFORMATION(1,equations_idx,variable_idx) = The equations type, see SOLVER_MAPPING_EquationsTypes
            !SUB_MATRIX_INFORMATION(2,equations_idx,variable_idx) = The equations set or interface condition index
            !SUB_MATRIX_INFORMATION(3,equations_idx,variable_idx) = The interface matrix index, or 0 for an equations set matrix
            !equations_idx goes from 1 to the number of equations sets + interface conditions
            !variable_idx goes from 1 to the number of variables mapped to this solver matrix
             ALLOCATE(subMatrixInformation(3,solverMapping%numberOfEquationsSets+solverMapping%numberOfInterfaceConditions, &
              & solverMapping%columnVariablesList(solverMatrixIdx)%numberOfVariables),STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate sub matrix information.",err,error,*999)
            subMatrixInformation=0
            !Allocate sub-matrix list information
            ALLOCATE(subMatrixList(0:3,solverMapping%numberOfEquationsSets+solverMapping%numberOfInterfaceConditions, &
              & solverMapping%columnVariablesList(solverMatrixIdx)%numberOfVariables),STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate sub matrix list.",err,error,*999)
            subMatrixList=0

            !
            ! 4b Calculate the number of columns
            !
            
            !Calculate the number of solver dofs
            ALLOCATE(numberOfVariableGlobalSolverDofs(solverMapping%columnVariablesList(solverMatrixIdx)%numberOfVariables), &
              & STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate number of global solver dofs.",err,error,*999)
            ALLOCATE(numberOfVariableLocalSolverDofs(solverMapping%columnVariablesList(solverMatrixIdx)%numberOfVariables), &
              & STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate number of local solver dofs.",err,error,*999)
            ALLOCATE(totalNumberOfVariableLocalSolverDofs(solverMapping%columnVariablesList(solverMatrixIdx)%numberOfVariables), &
              STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate total number of local solver dofs.",err,error,*999)
            
            numberOfVariableGlobalSolverDofs=0
            numberOfVariableLocalSolverDofs=0
            totalNumberOfVariableLocalSolverDofs=0
            
            equationsIdx=0
            !Loop over the equations sets
            DO equationsSetIdx=1,solverMapping%numberOfEquationsSets
              equationsIdx=equationsIdx+1
              !The pointers below have been checked for association above.
              equationsSet=>solverMapping%equationsSets(equationsSetIdx)%ptr
              equations=>equationsSet%equations
              equationsMapping=>equations%EQUATIONS_MAPPING
              dynamicMapping=>equationsMapping%DYNAMIC_MAPPING
              linearMapping=>equationsMapping%LINEAR_MAPPING
              nonlinearMapping=>equationsMapping%NONLINEAR_MAPPING
              dependentField=>equationsSet%dependent%DEPENDENT_FIELD
              boundaryConditions=>solverEquations%BOUNDARY_CONDITIONS
              NULLIFY(equationsSetVariableList)
              CALL List_CreateStart(equationsSetVariableList,err,error,*999)
              CALL List_DataTypeSet(equationsSetVariableList,LIST_INTG_TYPE,err,error,*999)
              CALL List_DataDimensionSet(equationsSetVariableList,3,err,error,*999)
              CALL List_KeyDimensionSet(equationsSetVariableList,1,err,error,*999)
              CALL List_CreateFinish(equationsSetVariableList,err,error,*999)
              IF(ASSOCIATED(dynamicMapping)) THEN
                equationsVariableListItem(1)=solverMapping%createValuesCache%dynamicVariableType(equationsSetIdx)
                equationsVariableListItem(2)=SOLVER_MAPPING_EQUATIONS_DYNAMIC_MATRIX
                equationsVariableListItem(3)=0
                CALL List_ItemAdd(equationsSetVariableList,equationsVariableListItem,err,error,*999)
              ENDIF
              IF(ASSOCIATED(linearMapping)) THEN
                DO variableIdx=1,solverMapping%createValuesCache%matrixVariableTypes(0,equationsSetIdx,solverMatrixIdx)
                  equationsVariableListItem(1)=solverMapping%createValuesCache%matrixVariableTypes( &
                    & variableIdx,equationsSetIdx,solverMatrixIdx)
                  equationsVariableListItem(2)=SOLVER_MAPPING_EQUATIONS_LINEAR_MATRIX
                  equationsVariableListItem(3)=variableIdx                
                  CALL List_ItemAdd(equationsSetVariableList,equationsVariableListItem,err,error,*999)
                ENDDO
              ENDIF
              IF(ASSOCIATED(nonlinearMapping)) THEN
                DO variableIdx=1,solverMapping%createValuesCache%residualVariableTypes(0,equationsSetIdx)
                  equationsVariableListItem(1)=solverMapping%createValuesCache%residualVariableTypes( &
                      & variableIdx,equationsSetIdx)
                  equationsVariableListItem(2)=SOLVER_MAPPING_EQUATIONS_NONLINEAR_MATRIX
                  equationsVariableListItem(3)=variableIdx
                  CALL List_ItemAdd(equationsSetVariableList,equationsVariableListItem,err,error,*999)
                ENDDO
              ENDIF
              CALL List_RemoveDuplicates(equationsSetVariableList,err,error,*999)
              CALL List_DetachAndDestroy(equationsSetVariableList,numberOfVariables,equationsSetVariables,err,error,*999)
              !Initialise equations set to solver map (sm)
              CALL SolverMapping_EquationsToSolverMatrixMapsSmInitialise(solverMapping% &
                & equationsSetToSolverMap(equationsSetIdx)%equationsToSolverMatrixMapsSm(solverMatrixIdx), &
                & err,error,*999)
              solverMapping%equationsSetToSolverMap(equationsSetIdx)%equationsToSolverMatrixMapsSm( &
                & solverMatrixIdx)%solverMatrixNumber=solverMatrixIdx
              !Allocate the equations set to solver map variables arrays
              ALLOCATE(solverMapping%equationsSetToSolverMap(equationsSetIdx)%equationsToSolverMatrixMapsSm( &
                & solverMatrixIdx)%variableTypes(numberOfVariables),STAT=err)
              IF(err/=0)  &
                & CALL FlagError("Could not allocate equations to solver matrix maps sm variable types.",err,error,*999)
              ALLOCATE(solverMapping%equationsSetToSolverMap(equationsSetIdx)%equationsToSolverMatrixMapsSm( &
                & solverMatrixIdx)%variables(numberOfVariables),STAT=err)
              IF(err/=0) CALL FlagError("Could not allocate equations to solver matrix maps sm variables.",err,error,*999)
              ALLOCATE(solverMapping%equationsSetToSolverMap(equationsSetIdx)%equationsToSolverMatrixMapsSm( &
                & solverMatrixIdx)%variableToSolverColMaps(numberOfVariables),STAT=err)
              IF(err/=0)  &
                & CALL FlagError("Could not allocate equations to solver matrix maps sm variables to solver col maps.", &
                & err,error,*999)                  
              !Setup
              solverMapping%equationsSetToSolverMap(equationsSetIdx)%equationsToSolverMatrixMapsSm( &
                & solverMatrixIdx)%numberOfVariables=numberOfVariables
              numberOfDynamicEquationsMatrices=0
              numberOfLinearEquationsMatrices=0
              !Loop over the variables in this equations set.
              DO variableIdx=1,numberOfVariables
                variableType=equationsSetVariables(1,variableIdx)
                matricesType=equationsSetVariables(2,variableIdx)
                matrixVariableIdx=equationsSetVariables(3,variableIdx)
                dependentVariable=>dependentField%VARIABLE_TYPE_MAP(variableType)%ptr
                IF(ASSOCIATED(dependentVariable)) THEN
                  !Find the variable in the list of solver variables
                  found=.FALSE.
                  DO variablePositionIdx=1,solverMapping%columnVariablesList(solverMatrixIdx)%numberOfVariables
                    IF(ASSOCIATED(dependentVariable,solverMapping%columnVariablesList(solverMatrixIdx)%variables( &
                      & variablePositionIdx)%variable)) THEN
                      found=.TRUE.
                      EXIT
                    ENDIF
                  ENDDO !variablePositionIdx
                  IF(found) THEN
                    !Add the equations set variable to the list of equations involving the solver variable
                    variableListItem(1)=equationsSetIdx
                    variableListItem(2)=variableType
                    variableListItem(3)=SOLVER_MAPPING_EQUATIONS_EQUATIONS_SET
                    CALL List_ItemAdd(variablesList(variablePositionIdx)%ptr,variableListItem,err,error,*999)
                    colDofsMapping=>dependentVariable%DOMAIN_MAPPING
                    IF(ASSOCIATED(colDofsMapping)) THEN
                      CALL BOUNDARY_CONDITIONS_VARIABLE_GET(boundaryConditions,dependentVariable,boundaryConditionsVariable, &
                        & err,error,*999)
                      IF(ASSOCIATED(boundaryConditionsVariable)) THEN
                        solverMapping%equationsSetToSolverMap(equationsSetIdx)%equationsToSolverMatrixMapsSm( &
                          & solverMatrixIdx)%variableTypes(variableIdx)=variableType
                        solverMapping%equationsSetToSolverMap(equationsSetIdx)%equationsToSolverMatrixMapsSm( &
                          & solverMatrixIdx)%variables(variableIdx)%ptr=>dependentVariable
                        !Allocate the variable to solver col maps arrays
                        ALLOCATE(solverMapping%equationsSetToSolverMap(equationsSetIdx)% &
                          & equationsToSolverMatrixMapsSm(solverMatrixIdx)%variableToSolverColMaps(variableIdx)% &
                          & columnNumbers(dependentVariable%TOTAL_NUMBER_OF_DOFS),STAT=err)
                        IF(err/=0) CALL FlagError("Could not allocate variables to solver column maps column numbers.", &
                          & err,error,*999)
                        ALLOCATE(solverMapping%equationsSetToSolverMap(equationsSetIdx)% &
                          & equationsToSolverMatrixMapsSm(solverMatrixIdx)%variableToSolverColMaps(variableIdx)% &
                          & couplingCoefficients(dependentVariable%TOTAL_NUMBER_OF_DOFS),STAT=err)
                        IF(err/=0) CALL FlagError("Could not allocate variables to solver column maps coupling coefficients.", &
                          & err,error,*999)
                        ALLOCATE(solverMapping%equationsSetToSolverMap(equationsSetIdx)% &
                          & equationsToSolverMatrixMapsSm(solverMatrixIdx)%variableToSolverColMaps(variableIdx)% &
                          & additiveConstants(dependentVariable%TOTAL_NUMBER_OF_DOFS),STAT=err)
                        IF(err/=0) CALL FlagError("Could not allocate variables to solver column maps additive constants.", &
                          & err,error,*999)
                        !Setup
                        !Set the sub-matrix information
                        subMatrixInformation(1,equationsIdx,variablePositionIdx)=SOLVER_MAPPING_EQUATIONS_EQUATIONS_SET
                        subMatrixInformation(2,equationsIdx,variablePositionIdx)=equationsSetIdx
                        !Set the sub-matrix lists
                        IF(ASSOCIATED(dynamicMapping)) THEN
                          numberOfDynamicEquationsMatrices=numberOfDynamicEquationsMatrices+ &
                            & dynamicMapping%VAR_TO_EQUATIONS_MATRICES_MAPS(variableType)%NUMBER_OF_EQUATIONS_MATRICES
                          IF(dynamicMapping%VAR_TO_EQUATIONS_MATRICES_MAPS(variableType)%NUMBER_OF_EQUATIONS_MATRICES>0) THEN
                            subMatrixList(0,equationsIdx,variablePositionIdx)= subMatrixList(0,equationsIdx,variablePositionIdx)+1
                            subMatrixList(subMatrixList(0,equationsIdx,variablePositionIdx), &
                              & equationsIdx,variablePositionIdx)=SOLVER_MAPPING_EQUATIONS_DYNAMIC_MATRIX
                          ENDIF
                        ENDIF
                        IF(ASSOCIATED(linearMapping)) THEN
                          numberOfLinearEquationsMatrices=numberOfLinearEquationsMatrices+ &
                            & linearMapping%VAR_TO_EQUATIONS_MATRICES_MAPS(variableType)%NUMBER_OF_EQUATIONS_MATRICES
                          IF(linearMapping%VAR_TO_EQUATIONS_MATRICES_MAPS(variableType)%NUMBER_OF_EQUATIONS_MATRICES>0) THEN
                            subMatrixList(0,equationsIdx,variablePositionIdx)=subMatrixList(0,equationsIdx,variablePositionIdx)+1
                            subMatrixList(subMatrixList(0,equationsIdx,variablePositionIdx), &
                              & equationsIdx,variablePositionIdx)=SOLVER_MAPPING_EQUATIONS_LINEAR_MATRIX
                          ENDIF
                        ENDIF
                        IF(ASSOCIATED(nonlinearMapping)) THEN
                          DO equationsMatrixIdx=1,nonlinearMapping%NUMBER_OF_RESIDUAL_VARIABLES
                            IF(nonlinearMapping%VAR_TO_JACOBIAN_MAP(equationsMatrixIdx)%VARIABLE_TYPE==variableType) THEN
                              subMatrixList(0,equationsIdx,variablePositionIdx)=subMatrixList(0,equationsIdx,variablePositionIdx)+1
                              subMatrixList(subMatrixList(0,equationsIdx,variablePositionIdx), &
                                & equationsIdx,variablePositionIdx)=SOLVER_MAPPING_EQUATIONS_NONLINEAR_MATRIX
                            ENDIF
                          ENDDO
                        ENDIF
                        !Loop over the global dofs for this variable.
                        DO globalDof=1,dependentVariable%NUMBER_OF_GLOBAL_DOFS
                          DO rankIdx=1,colDofsMapping%GLOBAL_TO_LOCAL_MAP(globalDof)%NUMBER_OF_DOMAINS
                            localDof=colDofsMapping%GLOBAL_TO_LOCAL_MAP(globalDof)%LOCAL_NUMBER(rankIdx)
                            dofType=colDofsMapping%GLOBAL_TO_LOCAL_MAP(globalDof)%LOCAL_TYPE(rankIdx)
                            columnRank=colDofsMapping%GLOBAL_TO_LOCAL_MAP(globalDof)%DOMAIN_NUMBER(rankIdx)
                            includeColumn=boundaryConditionsVariable%DOF_TYPES(globalDof)==BOUNDARY_CONDITION_DOF_FREE
                            constrainedDof=boundaryConditionsVariable%DOF_TYPES(globalDof)==BOUNDARY_CONDITION_DOF_CONSTRAINED
                            globalDofCouplingNumber=0
                            dofConstraints=>boundaryConditionsVariable%dofConstraints
                            IF(ASSOCIATED(dofConstraints)) THEN
                              IF(dofConstraints%numberOfConstraints>0) THEN
                                IF(ALLOCATED(dofConstraints%dofCouplings)) THEN
                                  IF(ASSOCIATED(dofConstraints%dofCouplings(globalDof)%ptr)) THEN
                                    CALL SolverDofCouplings_AddCoupling(columnCouplings, &
                                      & dofConstraints%dofCouplings(globalDof)%ptr, &
                                      & globalDofCouplingNumber,err,error,*999)
                                  END IF
                                ELSE
                                  CALL FlagError("DOF constraints DOF couplings are not allocated.",err,error,*999)
                                END IF
                              END IF
                            END IF
                            columnListItem(1)=globalDof
                            columnListItem(2)=localDof
                            columnListItem(5)=globalDofCouplingNumber
                            IF(dofType/=DOMAIN_LOCAL_GHOST) THEN
                              !DOF is not a ghost dof
                              IF(includeColumn) THEN
                                columnListItem(3)=1
                                IF(.NOT.variableProcessed(variablePositionIdx)) THEN
                                  numberOfVariableGlobalSolverDofs(variablePositionIdx)= &
                                    & numberOfVariableGlobalSolverDofs(variablePositionIdx)+1
                                  IF(columnRank==myrank) THEN
                                    numberOfVariableLocalSolverDofs(variablePositionIdx)= &
                                      & numberOfVariableLocalSolverDofs(variablePositionIdx)+1
                                    totalNumberOfVariableLocalSolverDofs(variablePositionIdx)= &
                                      & totalNumberOfVariableLocalSolverDofs(variablePositionIdx)+1
                                  ENDIF
                                ENDIF
                              ELSE IF(constrainedDof) THEN
                                columnListItem(3)=2
                              ELSE
                                columnListItem(3)=0
                              ENDIF
                              columnListItem(4)=variableIdx
                              CALL List_ItemAdd(rankGlobalColsLists(1,equationsIdx,variablePositionIdx,columnRank)%ptr, &
                                & columnListItem,err,error,*999)
                            ELSE
                              !DOF is a ghost dof
                              IF(includeColumn) THEN
                                columnListItem(3)=1
                                IF(.NOT.variableProcessed(variablePositionIdx)) THEN
                                  IF(columnRank==myrank) totalNumberOfVariableLocalSolverDofs(variablePositionIdx)= &
                                    & totalNumberOfVariableLocalSolverDofs(variablePositionIdx)+1
                                ENDIF
                              ELSE IF(constrainedDof) THEN
                                columnListItem(3)=2
                              ELSE
                                columnListItem(3)=0
                              ENDIF
                              columnListItem(4)=variableIdx
                              CALL List_ItemAdd(rankGlobalColsLists(2,equationsIdx,variablePositionIdx,columnRank)%ptr, &
                                & columnListItem,err,error,*999)
                            ENDIF
                          ENDDO !rankIdx
                        ENDDO !globalDof
                      ELSE
                        CALL FlagError("Boundary condition variable not associated.",err,error,*999)
                      ENDIF
                    ELSE
                      CALL FlagError("Equations matrix columns degree of freedom mapping is not associated.",err,error,*999)
                    ENDIF
                    variableProcessed(variablePositionIdx)=.TRUE.
                  ELSE
                    CALL FlagError("Dependent variable does not exist in the list of solver variables.",err,error,*999)
                  ENDIF
                ELSE
                  CALL FlagError("Dependent variable is not associated.",err,error,*999)
                ENDIF
              ENDDO !variableIdx
              IF(ALLOCATED(equationsSetVariables)) DEALLOCATE(equationsSetVariables)
            ENDDO !equationsSetIdx
            DO interfaceConditionIdx=1,solverMapping%numberOfInterfaceConditions
              equationsIdx=equationsIdx+1
              !The pointers below have been checked for association above.
              interfaceCondition=>solverMapping%interfaceConditions(interfaceConditionIdx)%ptr
              SELECT CASE(interfaceCondition%method)
              CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD)
                interfaceEquations=>interfaceCondition%INTERFACE_EQUATIONS
                interfaceMapping=>interfaceEquations%INTERFACE_MAPPING
                interfaceDependent=>interfaceCondition%dependent
                !Initialise interface condition to solver map (sm)
                CALL SolverMapping_InterfaceToSolverMatrixMapsSmInitialise(solverMapping%interfaceConditionToSolverMap( &
                  & interfaceConditionIdx)%interfaceToSolverMatrixMapsSm(solverMatrixIdx),err,error,*999)
                solverMapping%interfaceConditionToSolverMap(interfaceConditionIdx)%interfaceToSolverMatrixMapsSm( &
                  & solverMatrixIdx)%solverMatrixNumber=solverMatrixIdx
                !Allocate the interface to solver map variables arrays
                ALLOCATE(solverMapping%interfaceConditionToSolverMap(interfaceConditionIdx)% &
                  & interfaceToSolverMatrixMapsSm(solverMatrixIdx)%dependentVariableTypes( &
                  & interfaceMapping%NUMBER_OF_INTERFACE_MATRICES),STAT=err)
                IF(err/=0)  &
                  & CALL FlagError("Could not allocate interface to solver matrix maps sm dependent variable types.", &
                  & err,error,*999)
                ALLOCATE(solverMapping%interfaceConditionToSolverMap(interfaceConditionIdx)% &
                  & interfaceToSolverMatrixMapsSm(solverMatrixIdx)%dependentVariables(interfaceMapping% &
                  & NUMBER_OF_INTERFACE_MATRICES),STAT=err)
                IF(err/=0) CALL FlagError("Could not allocate equations to solver matrix maps sm variables.",err,error,*999)
                ALLOCATE(solverMapping%interfaceConditionToSolverMap(interfaceConditionIdx)% &
                  & interfaceToSolverMatrixMapsSm(solverMatrixIdx)%dependentVariableToSolverColMaps( &
                  & interfaceMapping%NUMBER_OF_INTERFACE_MATRICES),STAT=err)
                IF(err/=0) CALL FlagError("Could not allocate interface to solver matrix maps sm dependent variables "// &
                  & "to solver col maps.",err,error,*999)                  
                !First add in the Lagrange to solver variables
                lagrangeField=>interfaceCondition%lagrange%LAGRANGE_FIELD
                !\todo Lagrange variable type set to the first variable type for now
                variableType=1
                lagrangeVariable=>lagrangeField%VARIABLE_TYPE_MAP(variableType)%ptr
                IF(ASSOCIATED(lagrangeVariable)) THEN
                  !Setup
                  solverMapping%interfaceConditionToSolverMap(interfaceConditionIdx)%interfaceToSolverMatrixMapsSm( &
                    & solverMatrixIdx)%lagrangeVariableType=variableType
                  solverMapping%interfaceConditionToSolverMap(interfaceConditionIdx)%interfaceToSolverMatrixMapsSm( &
                    & solverMatrixIdx)%lagrangeVariable=>lagrangeField%VARIABLE_TYPE_MAP(variableType)%ptr
                  !Find the variable in the list of solver variables
                  found=.FALSE.
                  DO variablePositionIdx=1,solverMapping%columnVariablesList(solverMatrixIdx)%numberOfVariables
                    IF(ASSOCIATED(lagrangeVariable,solverMapping%columnVariablesList(solverMatrixIdx)%variables( &
                      & variablePositionIdx)%variable)) THEN
                      found=.TRUE.
                      EXIT
                    ENDIF
                  ENDDO !variablePositionIdx
                  IF(found) THEN
                    !Add the interface condition variable to the list of equations involving the solver variable
                    variableListItem(1)=interfaceConditionIdx
                    variableListItem(2)=variableType
                    variableListItem(3)=SOLVER_MAPPING_EQUATIONS_INTERFACE_CONDITION
                    CALL List_ItemAdd(variablesList(variablePositionIdx)%ptr,variableListItem,err,error,*999)
                    colDofsMapping=>lagrangeVariable%DOMAIN_MAPPING
                    IF(ASSOCIATED(colDofsMapping)) THEN
                      CALL BOUNDARY_CONDITIONS_VARIABLE_GET(boundaryConditions,lagrangeVariable,boundaryConditionsVariable, &
                        & err,error,*999)
                      IF(ASSOCIATED(boundaryConditionsVariable)) THEN
                        solverMapping%interfaceConditionToSolverMap(interfaceConditionIdx)% &
                          & interfaceToSolverMatrixMapsSm(solverMatrixIdx)%lagrangeVariableType=variableType
                        solverMapping%interfaceConditionToSolverMap(interfaceConditionIdx)% &
                          & interfaceToSolverMatrixMapsSm(solverMatrixIdx)%lagrangeVariable=>lagrangeVariable
                        !Allocate the variable to solver col maps arrays
                        ALLOCATE(solverMapping%interfaceConditionToSolverMap(interfaceConditionIdx)% &
                          & interfaceToSolverMatrixMapsSm(solverMatrixIdx)%lagrangeVariableToSolverColMap% &
                          & columnNumbers(lagrangeVariable%TOTAL_NUMBER_OF_DOFS),STAT=err)
                        IF(err/=0) CALL FlagError("Could not allocate variables to solver column maps column numbers.", &
                          & err,error,*999)
                        ALLOCATE(solverMapping%interfaceConditionToSolverMap(interfaceConditionIdx)% &
                          & interfaceToSolverMatrixMapsSm(solverMatrixIdx)%lagrangeVariableToSolverColMap% &
                          & couplingCoefficients(lagrangeVariable%TOTAL_NUMBER_OF_DOFS),STAT=err)
                        IF(err/=0) CALL FlagError("Could not allocate variables to solver column maps coupling coefficients.", &
                          & err,error,*999)
                        ALLOCATE(solverMapping%interfaceConditionToSolverMap(interfaceConditionIdx)% &
                        & interfaceToSolverMatrixMapsSm(solverMatrixIdx)%lagrangeVariableToSolverColMap% &
                        & additiveConstants(lagrangeVariable%TOTAL_NUMBER_OF_DOFS),STAT=err)
                        IF(err/=0) CALL FlagError("Could not allocate variables to solver column maps additive constants.", &
                          & err,error,*999)
                        DO equationsIdx2=1,solverMapping%interfaceConditionToSolverMap(interfaceConditionIdx)% &
                          & numberOfEquationsSets
                          equationsSetIdx=solverMapping%interfaceConditionToSolverMap(interfaceConditionIdx)% &
                            & interfaceToSolverMatrixMapsEquations(equationsIdx2)%equationsSetIndex
                          interfaceMatrixIdx=solverMapping%interfaceConditionToSolverMap(interfaceConditionIdx)% &
                            & interfaceToSolverMatrixMapsEquations(equationsIdx2)%interfaceMatrixIndex
                          !Set the sub-matrix information
                          subMatrixInformation(1,equationsSetIdx,variablePositionIdx)=SOLVER_MAPPING_EQUATIONS_INTERFACE_CONDITION
                          subMatrixInformation(2,equationsSetIdx,variablePositionIdx)=interfaceConditionIdx
                          subMatrixInformation(3,equationsSetIdx,variablePositionIdx)=interfaceMatrixIdx
                          !Loop over the global dofs for this variable.
                          DO globalDof=1,lagrangeVariable%NUMBER_OF_GLOBAL_DOFS
                            DO rankIdx=1,colDofsMapping%GLOBAL_TO_LOCAL_MAP(globalDof)%NUMBER_OF_DOMAINS
                              localDof=colDofsMapping%GLOBAL_TO_LOCAL_MAP(globalDof)%LOCAL_NUMBER(rankIdx)
                              dofType=colDofsMapping%GLOBAL_TO_LOCAL_MAP(globalDof)%LOCAL_TYPE(rankIdx)
                              columnRank=colDofsMapping%GLOBAL_TO_LOCAL_MAP(globalDof)%DOMAIN_NUMBER(rankIdx)
                              !For now include Lagrange column
                              includeColumn=.TRUE.
                              columnListItem(1)=globalDof
                              columnListItem(2)=localDof
                              columnListItem(5)=0
                              IF(dofType/=DOMAIN_LOCAL_GHOST) THEN
                                !DOF is not a ghost dof
                                IF(includeColumn) THEN
                                  columnListItem(3)=1
                                  IF(.NOT.variableProcessed(variablePositionIdx)) THEN
                                    numberOfVariableGlobalSolverDofs(variablePositionIdx)= &
                                      & numberOfVariableGlobalSolverDofs(variablePositionIdx)+1
                                    IF(columnRank==myrank) THEN
                                      numberOfVariableLocalSolverDofs(variablePositionIdx)= &
                                        & numberOfVariableLocalSolverDofs(variablePositionIdx)+1
                                      totalNumberOfVariableLocalSolverDofs(variablePositionIdx)= &
                                        & totalNumberOfVariableLocalSolverDofs(variablePositionIdx)+1
                                    ENDIF
                                  ELSE
                                    columnListItem(3)=0
                                  ENDIF
                                  columnListItem(4)=variableIdx
                                  CALL List_ItemAdd(rankGlobalColsLists(1,equationsSetIdx,variablePositionIdx, &
                                    & columnRank)%ptr,columnListItem,err,error,*999)
                                ELSE
                                  columnListItem(3)=0
                                ENDIF
                                columnListItem(4)=variableIdx
                                CALL List_ItemAdd(rankGlobalColsLists(1,equationsSetIdx,variablePositionIdx, &
                                  & columnRank)%ptr,columnListItem,err,error,*999)
                              ELSE
                                !DOF is a ghost dof
                                IF(includeColumn) THEN
                                  columnListItem(3)=1
                                  IF(.NOT.variableProcessed(variablePositionIdx)) THEN
                                    IF(columnRank==myrank) totalNumberOfVariableLocalSolverDofs(variablePositionIdx)= &
                                      & totalNumberOfVariableLocalSolverDofs(variablePositionIdx)+1
                                  ENDIF
                                ELSE
                                  columnListItem(3)=0
                                ENDIF
                                columnListItem(4)=variableIdx
                                CALL List_ItemAdd(rankGlobalColsLists(2,equationsSetIdx,variablePositionIdx, &
                                  & columnRank)%ptr,columnListItem,err,error,*999)
                              ENDIF
                            ENDDO !rankIdx
                          ENDDO !globalDof
                          variableProcessed(variablePositionIdx)=.TRUE.
                        ENDDO !equationsIdx2
                      ELSE
                        CALL FLAG_ERROR("Boundary condition variable not associated.",ERR,ERROR,*999)
                      ENDIF
                    ELSE
                      CALL FlagError("Columns degree of freedom mapping is not associated.",err,error,*999)
                    ENDIF
                  ELSE
                    CALL FlagError("Lagrange variable does not exist in the list of solver variables.",err,error,*999)
                  ENDIF
                ELSE
                  CALL FlagError("Lagrange variable is not associated.",err,error,*999)
                ENDIF
                !Now add in the Dependent variables
                solverMapping%interfaceConditionToSolverMap(interfaceConditionIdx)%interfaceToSolverMatrixMapsSm( &
                  & solverMatrixIdx)%numberOfDependentVariables=interfaceMapping%NUMBER_OF_INTERFACE_MATRICES
                DO interfaceMatrixIdx=1,interfaceMapping%NUMBER_OF_INTERFACE_MATRICES
                  dependentVariable=>interfaceMapping%INTERFACE_MATRIX_ROWS_TO_VAR_MAPS(interfaceMatrixIdx)%variable
                  IF(ASSOCIATED(dependentVariable)) THEN
                    variableType=dependentVariable%VARIABLE_TYPE
                    !Find the variable in the list of solver variables
                    found=.FALSE.
                    DO variablePositionIdx=1,solverMapping%columnVariablesList(solverMatrixIdx)%numberOfVariables
                      IF(ASSOCIATED(dependentVariable,solverMapping%columnVariablesList(solverMatrixIdx)%variables( &
                        & variablePositionIdx)%VARIABLE)) THEN
                        found=.TRUE.
                        EXIT
                      ENDIF
                    ENDDO !variablePositionIdx
                    IF(found) THEN
                      equationsSet=>interfaceDependent%EQUATIONS_SETS(interfaceMapping%INTERFACE_MATRIX_ROWS_TO_VAR_MAPS( &
                        & interfaceMatrixIdx)%MESH_INDEX)%ptr
                      !Note that EQUATIONS_SET and INTERFACE_EQUATIONS has already been checked for association above and this check is just to see if either an equation set or interface equations are present.
                      IF(ASSOCIATED(equationsSet).OR.ASSOCIATED(interfaceEquations)) THEN
                          colDofsMapping=>dependentVariable%DOMAIN_MAPPING
                          IF(ASSOCIATED(colDofsMapping)) THEN
                            CALL BOUNDARY_CONDITIONS_VARIABLE_GET(boundaryConditions,dependentVariable, &
                              & boundaryConditionsVariable,err,error,*999)
                            IF(ASSOCIATED(boundaryConditionsVariable)) THEN
                              !Setup
                              solverMapping%interfaceConditionToSolverMap(interfaceConditionIdx)% &
                                & interfaceToSolverMatrixMapsSm(solverMatrixIdx)%dependentVariableTypes( &
                                & interfaceMatrixIdx)=variableType
                              solverMapping%interfaceConditionToSolverMap(interfaceConditionIdx)% &
                                & interfaceToSolverMatrixMapsSm(solverMatrixIdx)%dependentVariables(interfaceMatrixIdx)% &
                                & ptr=>dependentVariable
                              ALLOCATE(solverMapping%interfaceConditionToSolverMap(interfaceConditionIdx)% &
                                & interfaceToSolverMatrixMapsSm(solverMatrixIdx)%dependentVariableToSolverColMaps( &
                                & interfaceMatrixIdx)%columnNumbers(dependentVariable%TOTAL_NUMBER_OF_DOFS),STAT=err)
                              IF(err/=0) CALL FlagError("Could not allocate variables to solver column maps column numbers.", &
                                & err,error,*999)
                              ALLOCATE(solverMapping%interfaceConditionToSolverMap(interfaceConditionIdx)% &
                                & interfaceToSolverMatrixMapsSm(solverMatrixIdx)%dependentVariableToSolverColMaps( &
                                & interfaceMatrixIdx)%couplingCoefficients(dependentVariable%TOTAL_NUMBER_OF_DOFS),STAT=err)
                              IF(err/=0)  &
                                & CALL FlagError("Could not allocate variables to solver column maps coupling coefficients.", &
                                & err,error,*999)
                              ALLOCATE(solverMapping%interfaceConditionToSolverMap(interfaceConditionIdx)% &
                                & interfaceToSolverMatrixMapsSm(solverMatrixIdx)%dependentVariableToSolverColMaps( &
                                & interfaceMatrixIdx)%additiveConstants(dependentVariable%TOTAL_NUMBER_OF_DOFS),STAT=err)
                              IF(err/=0) CALL FlagError("Could not allocate variables to solver column maps additive constants.", &
                                & err,error,*999)
                              !Set the sub-matrix information
                              subMatrixInformation(1,equationsIdx,variablePositionIdx)=SOLVER_MAPPING_EQUATIONS_INTERFACE_TRANSPOSE
                              subMatrixInformation(2,equationsIdx,variablePositionIdx)=interfaceConditionIdx
                              subMatrixInformation(3,equationsIdx,variablePositionIdx)=interfaceMatrixIdx
                              !Loop over the global dofs for this variable.
                              DO globalDof=1,dependentVariable%NUMBER_OF_GLOBAL_DOFS
                                DO rankIdx=1,colDofsMapping%GLOBAL_TO_LOCAL_MAP(globalDof)%NUMBER_OF_DOMAINS
                                  localDof=colDofsMapping%GLOBAL_TO_LOCAL_MAP(globalDof)%LOCAL_NUMBER(rankIdx)
                                  dofType=colDofsMapping%GLOBAL_TO_LOCAL_MAP(globalDof)%LOCAL_TYPE(rankIdx)
                                  columnRank=colDofsMapping%GLOBAL_TO_LOCAL_MAP(globalDof)%DOMAIN_NUMBER(rankIdx)
                                  includeColumn=boundaryConditionsVariable%DOF_TYPES(globalDof)==BOUNDARY_CONDITION_DOF_FREE
                                  constrainedDof=boundaryConditionsVariable%DOF_TYPES(globalDof)==BOUNDARY_CONDITION_DOF_CONSTRAINED
                                  globalDofCouplingNumber=0
                                  dofConstraints=>boundaryConditionsVariable%dofConstraints
                                  IF(ASSOCIATED(dofConstraints)) THEN
                                    IF(dofConstraints%numberOfConstraints>0) THEN
                                      IF(ALLOCATED(dofConstraints%dofCouplings)) THEN
                                        IF(ASSOCIATED(dofConstraints%dofCouplings(globalDof)%ptr)) THEN
                                          CALL SolverDofCouplings_AddCoupling(columnCouplings, &
                                            & dofConstraints%dofCouplings(globalDof)%ptr, &
                                            & globalDofCouplingNumber,err,error,*999)
                                        END IF
                                      ELSE
                                        CALL FlagError("DOF constraints DOF couplings are not allocated.",err,error,*999)
                                      END IF
                                    END IF
                                  END IF
                                  columnListItem(1)=globalDof
                                  columnListItem(2)=localDof
                                  columnListItem(5)=globalDofCouplingNumber
                                  IF(dofType/=DOMAIN_LOCAL_GHOST) THEN
                                    !DOF is not a ghost dof
                                    IF(includeColumn) THEN
                                      columnListItem(3)=1
                                      IF(.NOT.variableProcessed(variablePositionIdx)) THEN
                                        numberOfVariableGlobalSolverDofs(variablePositionIdx)= &
                                          & numberOfVariableGlobalSolverDofs(variablePositionIdx)+1
                                        IF(columnRank==myrank) THEN
                                          numberOfVariableLocalSolverDofs(variablePositionIdx)= &
                                            & numberOfVariableLocalSolverDofs(variablePositionIdx)+1
                                          totalNumberOfVariableLocalSolverDofs(variablePositionIdx)= &
                                            & totalNumberOfVariableLocalSolverDofs(variablePositionIdx)+1
                                        ENDIF
                                      ENDIF
                                    ELSE IF(constrainedDof) THEN
                                      columnListItem(3)=2
                                    ELSE
                                      columnListItem(3)=0
                                    ENDIF
                                    columnListItem(4)=variableIdx
                                    CALL List_ItemAdd(rankGlobalColsLists(1,equationsIdx,variablePositionIdx,columnRank)%ptr, &
                                      & columnListItem,err,error,*999)
                                  ELSE
                                    !DOF is a ghost dof
                                    IF(includeColumn) THEN
                                      columnListItem(3)=1
                                      IF(.NOT.variableProcessed(variablePositionIdx)) THEN
                                        IF(columnRank==myrank) totalNumberOfVariableLocalSolverDofs(variablePositionIdx)= &
                                          & totalNumberOfVariableLocalSolverDofs(variablePositionIdx)+1
                                      ENDIF
                                    ELSE IF(constrainedDof) THEN
                                      columnListItem(3)=2
                                    ELSE
                                      columnListItem(3)=0
                                    ENDIF
                                    columnListItem(4)=variableIdx
                                    CALL List_ItemAdd(rankGlobalColsLists(2,equationsIdx,variablePositionIdx, &
                                      & columnRank)%ptr,columnListItem,err,error,*999)
                                  ENDIF
                                ENDDO !rankIdx
                              ENDDO !globalDof
                            ELSE
                              CALL FlagError("Boundary condition variable not associated.",err,error,*999)
                            ENDIF
                          ELSE
                            CALL FlagError("Interface matrix columns degree of freedom mapping is not associated.",err,error,*999)
                          ENDIF
                        ELSE
                          CALL FlagError("Interface depdendent equations set is not associated.",err,error,*999)
                      ENDIF
                      variableProcessed(variablePositionIdx)=.TRUE.
                    ELSE
                      CALL FlagError("Dependent variable does not exist in the list of solver variables.",err,error,*999)
                    ENDIF
                  ELSE
                    CALL FlagError("Dependent variable is not associated.",err,error,*999)
                  ENDIF
                ENDDO !matrix_idx
              CASE(INTERFACE_CONDITION_AUGMENTED_LAGRANGE_METHOD)
                CALL FlagError("Not implemented.",err,error,*999)
              CASE(INTERFACE_CONDITION_PENALTY_METHOD)
                CALL FlagError("Not implemented.",err,error,*999)
              CASE(INTERFACE_CONDITION_POINT_TO_POINT_METHOD)
                CALL FlagError("Not implemented.",err,error,*999)
              CASE DEFAULT
                localError="The interface condition method of "// &
                  & TRIM(NumberToVString(interfaceCondition%method,"*",err,error))// &
                  & " is invalid."
                CALL FlagError(localError,err,error,*999)
              END SELECT
            ENDDO !interfaceIdx

            IF(ALLOCATED(variableProcessed)) DEALLOCATE(variableProcessed)

            numberOfLocalSolverDofs=0
            totalNumberOfLocalSolverDofs=0
            numberOfGlobalSolverDofs=0
            DO solverVariableIdx=1,solverMapping%columnVariablesList(solverMatrixIdx)%numberOfVariables
              numberOfLocalSolverDofs=numberOfLocalSolverDofs+numberOfVariableLocalSolverDofs(solverVariableIdx)
              totalNumberOfLocalSolverDofs=totalNumberOfLocalSolverDofs+totalNumberOfVariableLocalSolverDofs(solverVariableIdx)
              numberOfGlobalSolverDofs=numberOfGlobalSolverDofs+numberOfVariableGlobalSolverDofs(solverVariableIdx)
            ENDDO !solverVariableIdx

            !Sanity check
            IF(numberOfLocalSolverDofs==0) THEN
              localError="Invalid problem setup. The number of local solver DOFs for solver matrix "// &
                & TRIM(NumberToVString(solverMatrixIdx,"*",err,error))//" is zero."
              CALL FlagError(localError,err,error,*999)
            ENDIF
            IF(numberOfGlobalSolverDofs==0) THEN
              localError="Invalid problem setup. The number of global solver DOFs for solver matrix "// &
                & TRIM(NumberToVString(solverMatrixIdx,"*",err,error))//" is zero."
              CALL FlagError(localError,err,error,*999)
            ENDIF

            !
            ! 4c Set up, allocate and initialise column mappings
            !

            !Allocate memory for this solver matrix
            !Allocate solver columns to equations sets maps
            ALLOCATE(solverMapping%solverColToEquationsColsMap(solverMatrixIdx)%solverDofToVariableMaps( &
              & totalNumberOfLocalSolverDofs),STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate solver dof to variable maps.",err,error,*999)
            !Set the number of columns
            solverMapping%solverColToEquationsColsMap(solverMatrixIdx)%numberOfColumns=numberOfGlobalSolverDofs
            !Set the number of variables
            solverMapping%solverColToEquationsColsMap(solverMatrixIdx)%numberOfDofs=numberOfLocalSolverDofs
            solverMapping%solverColToEquationsColsMap(solverMatrixIdx)%totalNumberOfDofs=totalNumberOfLocalSolverDofs
            solverMapping%solverColToEquationsColsMap(solverMatrixIdx)%numberOfGlobalDofs=numberOfGlobalSolverDofs
            !Allocate the columns domain mapping
            ALLOCATE(solverMapping%solverColToEquationsColsMap(solverMatrixIdx)%columnDofsMapping,STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate solver col to equations sets map column dofs mapping.",err,error,*999)
!!TODO: what is the real number of domains for a solver???
            CALL DOMAIN_MAPPINGS_MAPPING_INITIALISE(solverMapping%solverColToEquationsColsMap(solverMatrixIdx)% &
              & columnDofsMapping,COMPUTATIONAL_ENVIRONMENT%NUMBER_COMPUTATIONAL_NODES,err,error,*999)            
            colDomainMapping=>solverMapping%solverColToEquationsColsMap(solverMatrixIdx)%columnDofsMapping
            ALLOCATE(colDomainMapping%GLOBAL_TO_LOCAL_MAP(numberOfGlobalSolverDofs),STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate column dofs mapping global to local.",err,error,*999)
            colDomainMapping%NUMBER_OF_GLOBAL=numberOfGlobalSolverDofs
            ALLOCATE(variableRankProcessed(solverMapping%columnVariablesList(solverMatrixIdx)%numberOfVariables, &
              & 0:COMPUTATIONAL_ENVIRONMENT%NUMBER_COMPUTATIONAL_NODES-1),STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate variable rank processed.",err,error,*999)
            variableRankProcessed=.FALSE.
            !Calculate the column mappings
            numberOfColumns=numberOfGlobalSolverDofs

            !Initialise            
            DO equationsSetIdx=1,solverMapping%numberOfEquationsSets

              equationsSet=>solverMapping%equationsSets(equationsSetIdx)%ptr
              equations=>equationsSet%equations
              equationsMapping=>equations%EQUATIONS_MAPPING
              dynamicMapping=>equationsMapping%DYNAMIC_MAPPING
              linearMapping=>equationsMapping%LINEAR_MAPPING
              nonlinearMapping=>equationsMapping%NONLINEAR_MAPPING
              
              IF(ASSOCIATED(dynamicMapping)) THEN
                !Allocate the equations set to solver maps for equations matrix (em) indexing
                ALLOCATE(solverMapping%equationsSetToSolverMap(equationsSetIdx)%equationsToSolverMatrixMapsEm( &
                  & dynamicMapping%NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES),STAT=err)
                IF(err/=0) &
                  & CALL FlagError("Could not allocate equations set to solver map equations to solver matrix maps em.", &
                  & err,error,*999)
                DO equationsMatrixIdx=1,dynamicMapping%NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES
                  CALL SolverMapping_EquationsToSolverMatrixMapsEmInitialise(solverMapping%equationsSetToSolverMap( &
                    & equationsSetIdx)%equationsToSolverMatrixMapsEm(equationsMatrixIdx),err,error,*999)
                ENDDO !equationsMatrixIdx
                IF(ASSOCIATED(nonlinearMapping)) THEN
                  !Allocate the equations set to solver maps for Jacobian matrix (jm) indexing
                  ALLOCATE(solverMapping%equationsSetToSolverMap(equationsSetIdx)%equationsToSolverMatrixMapsJm( &
                    & nonlinearMapping%NUMBER_OF_RESIDUAL_VARIABLES),STAT=err)
                  IF(err/=0) &
                    & CALL FlagError("Could not allocate equations set to solver map equations to solver matrix maps jm.", &
                    & err,error,*999)
                  DO equationsMatrixIdx=1,nonlinearMapping%NUMBER_OF_RESIDUAL_VARIABLES
                    NULLIFY(solverMapping%equationsSetToSolverMap(equationsSetIdx)%equationsToSolverMatrixMapsJm( &
                      & equationsMatrixIdx)%ptr)
                  ENDDO
                ENDIF
              ELSE
                IF(ASSOCIATED(linearMapping)) THEN
                  !Allocate the equations set to solver maps for equations matrix (em) indexing
                  ALLOCATE(solverMapping%equationsSetToSolverMap(equationsSetIdx)%equationsToSolverMatrixMapsEm( &
                    & linearMapping%NUMBER_OF_LINEAR_EQUATIONS_MATRICES),STAT=err)
                  IF(err/=0) &
                    & CALL FlagError("Could not allocate equations set to solver map equations to solver matrix maps em.", &
                    & err,error,*999)
                  DO equationsMatrixIdx=1,linearMapping%NUMBER_OF_LINEAR_EQUATIONS_MATRICES
                    CALL SolverMapping_EquationsToSolverMatrixMapsEmInitialise(solverMapping% &
                      & equationsSetToSolverMap(equationsSetIdx)%equationsToSolverMatrixMapsEm(equationsMatrixIdx), &
                      & err,error,*999)
                  ENDDO !equationsMatrixIdx
                ENDIF
                IF(ASSOCIATED(nonlinearMapping)) THEN
                  !Allocate the equations set to solver maps for Jacobian matrix (jm) indexing
                  ALLOCATE(solverMapping%equationsSetToSolverMap(equationsSetIdx)%equationsToSolverMatrixMapsJm( &
                    & nonlinearMapping%NUMBER_OF_RESIDUAL_VARIABLES),STAT=err)
                  IF(err/=0) &
                    & CALL FlagError("Could not allocate equations set to solver map equations to solver matrix maps jm.", &
                    & err,error,*999)
                  DO equationsMatrixIdx=1,nonlinearMapping%NUMBER_OF_RESIDUAL_VARIABLES
                    NULLIFY(solverMapping%equationsSetToSolverMap(equationsSetIdx)%equationsToSolverMatrixMapsJm( &
                      & equationsMatrixIdx)%ptr)
                  ENDDO
                ENDIF
              ENDIF
              
              !Initialise solver columns to equations set map
              CALL SolverMapping_SolverColToEquationsSetMapInitialise(solverMapping% &
                & solverColToEquationsColsMap(solverMatrixIdx)%solverColToEquationsSetMaps( &
                & equationsSetIdx),err,error,*999)
              
              dependentField=>equationsSet%dependent%DEPENDENT_FIELD
              IF(ASSOCIATED(dynamicMapping)) THEN
                numberOfVariables=1
              ELSE
                IF(ASSOCIATED(nonlinearMapping)) THEN
                  numberOfVariables=solverMapping%createValuesCache%residualVariableTypes(0,equationsSetIdx)
                ELSE
                  numberOfVariables=solverMapping%createValuesCache%matrixVariableTypes(0,equationsSetIdx,solverMatrixIdx)
                ENDIF
              ENDIF
              
              !Allocate the solver columns to equations set map arrays
              IF(ASSOCIATED(dynamicMapping)) THEN
                solverMapping%solverColToEquationsColsMap(solverMatrixIdx)%solverColToEquationsSetMaps( &
                  & equationsSetIdx)%haveDynamic=.TRUE.
                ALLOCATE(solverMapping%solverColToEquationsColsMap(solverMatrixIdx)%solverColToEquationsSetMaps( &
                  & equationsSetIdx)%solverColToDynamicEquationsMaps(numberOfColumns),STAT=err)
                IF(err/=0) CALL FlagError("Could not allocate solver columns to dynamic equations map.",err,error,*999)
              ELSE
                solverMapping%solverColToEquationsColsMap(solverMatrixIdx)%solverColToEquationsSetMaps( &
                  & equationsSetIdx)%haveStatic=.TRUE.
                ALLOCATE(solverMapping%solverColToEquationsColsMap(solverMatrixIdx)%solverColToEquationsSetMaps( &
                  & equationsSetIdx)%SolverColToStaticEquationsMaps(numberOfColumns),STAT=err)
                IF(err/=0) CALL FlagError("Could not allocate solver columns to static equations map.",err,error,*999)
              ENDIF
              !Set the solver column to equations set map
              solverMapping%solverColToEquationsColsMap(solverMatrixIdx)%solverColToEquationsSetMaps( &
                & equationsSetIdx)%equations=>equations
              
              !Allocate the equations to solver matrix maps sm equations to solver maps
              IF(ASSOCIATED(dynamicMapping)) THEN
                ALLOCATE(solverMapping%equationsSetToSolverMap(equationsSetIdx)%equationsToSolverMatrixMapsSm( &
                  & solverMatrixIdx)%dynamicEquationsToSolverMatrixMaps(numberOfDynamicEquationsMatrices),STAT=err)
                IF(err/=0) CALL FlagError("Could not allocate equations to solver matrix maps sm dynamic equations "// &
                  & "to solver matrix maps.",err,error,*999)
                !Set up dynamic arrays
                solverMapping%equationsSetToSolverMap(equationsSetIdx)%equationsToSolverMatrixMapsSm( &
                  & solverMatrixIdx)%numberOfDynamicEquationsMatrices=numberOfDynamicEquationsMatrices
                DO equationsMatrixIdx=1,numberOfDynamicEquationsMatrices
                  NULLIFY(solverMapping%equationsSetToSolverMap(equationsSetIdx)%equationsToSolverMatrixMapsSm( &
                    & solverMatrixIdx)%dynamicEquationsToSolverMatrixMaps(equationsMatrixIdx)%ptr)
                  ALLOCATE(solverMapping%equationsSetToSolverMap(equationsSetIdx)%equationsToSolverMatrixMapsSm( &
                    & solverMatrixIdx)%dynamicEquationsToSolverMatrixMaps(equationsMatrixIdx)%ptr,STAT=err)
                  IF(err/=0) CALL FlagError("Could not allocate equations to solver matrix maps.",err,error,*999)
                  CALL SolverMapping_EquationsToSolverMapsInitialise(solverMapping%equationsSetToSolverMap( &
                    & equationsSetIdx)%equationsToSolverMatrixMapsSm(solverMatrixIdx)% &
                    & dynamicEquationsToSolverMatrixMaps(equationsMatrixIdx)%ptr,err,error,*999)
                  solverMapping%equationsSetToSolverMap(equationsSetIdx)%equationsToSolverMatrixMapsSm( &
                    & solverMatrixIdx)%dynamicEquationsToSolverMatrixMaps(equationsMatrixIdx)%ptr% &
                    & equationsMatrixType=SOLVER_MAPPING_EQUATIONS_DYNAMIC_MATRIX
                  solverMapping%equationsSetToSolverMap(equationsSetIdx)%equationsToSolverMatrixMapsEm( &
                    & equationsMatrixIdx)%equationsMatrixNumber=equationsMatrixIdx
                  solverMapping%equationsSetToSolverMap(equationsSetIdx)%equationsToSolverMatrixMapsEm( &
                    & equationsMatrixIdx)%numberOfSolverMatrices=solverMapping%equationsSetToSolverMap( &
                    & equationsSetIdx)%equationsToSolverMatrixMapsEm(equationsMatrixIdx)%numberOfSolverMatrices+1
                ENDDO !equationsMatrixIdx                    
                !Set up nonlinear arrays
                IF(ASSOCIATED(nonlinearMapping)) THEN
                  solverMapping%equationsSetToSolverMap(equationsSetIdx)%equationsToSolverMatrixMapsSm( &
                    & solverMatrixIdx)%numberOfEquationsJacobians=nonlinearMapping%NUMBER_OF_RESIDUAL_VARIABLES
                  ALLOCATE(solverMapping%equationsSetToSolverMap(equationsSetIdx)%equationsToSolverMatrixMapsSm( &
                    & solverMatrixIdx)%jacobianToSolverMatrixMaps(nonlinearMapping%NUMBER_OF_RESIDUAL_VARIABLES),STAT=err)
                  IF(err/=0) CALL FlagError("Could not allocate Jacobian to solver matrix maps.",err,error,*999)
                  DO equationsMatrixIdx=1,nonlinearMapping%NUMBER_OF_RESIDUAL_VARIABLES
                    ALLOCATE(solverMapping%equationsSetToSolverMap(equationsSetIdx)%equationsToSolverMatrixMapsSm( &
                      & solverMatrixIdx)%jacobianToSolverMatrixMaps(equationsMatrixIdx)%ptr,STAT=err)
                    IF(err/=0) CALL FlagError("Could not allocate Jacobian to solver matrix map.",err,error,*999)
                    CALL SolverMapping_JacobianToSolverMapInitialise(solverMapping%equationsSetToSolverMap( &
                      & equationsSetIdx)%equationsToSolverMatrixMapsSm(solverMatrixIdx)% &
                      & jacobianToSolverMatrixMaps(equationsMatrixIdx)%ptr,err,error,*999)
                    solverMapping%equationsSetToSolverMap(equationsSetIdx)%equationsToSolverMatrixMapsJm &
                      & (equationsMatrixIdx)%ptr=>solverMapping%equationsSetToSolverMap(equationsSetIdx)% &
                      & equationsToSolverMatrixMapsSm(solverMatrixIdx)%jacobianToSolverMatrixMaps( &
                      & equationsMatrixIdx)%ptr
                  ENDDO !equationsMatrixIdx
                ENDIF
              ELSE
                IF(ASSOCIATED(linearMapping)) THEN
                  ALLOCATE(solverMapping%equationsSetToSolverMap(equationsSetIdx)%equationsToSolverMatrixMapsSm( &
                    & solverMatrixIdx)%linearEquationsToSolverMatrixMaps(numberOfLinearEquationsMatrices),STAT=err)
                  IF(err/=0) &
                    & CALL FlagError("Could not allocate equations to solver matrix maps sm equations to solver matrix maps.", &
                    & err,error,*999)
                  !Set up linear arrays
                  solverMapping%equationsSetToSolverMap(equationsSetIdx)%equationsToSolverMatrixMapsSm( &
                    & solverMatrixIdx)%numberOfLinearEquationsMatrices=numberOfLinearEquationsMatrices
                  DO equationsMatrixIdx=1,numberOfLinearEquationsMatrices
                    NULLIFY(solverMapping%equationsSetToSolverMap(equationsSetIdx)%equationsToSolverMatrixMapsSm( &
                      & solverMatrixIdx)%linearEquationsToSolverMatrixMaps(equationsMatrixIdx)%ptr)
                    ALLOCATE(solverMapping%equationsSetToSolverMap(equationsSetIdx)%equationsToSolverMatrixMapsSm( &
                      & solverMatrixIdx)%linearEquationsToSolverMatrixMaps(equationsMatrixIdx)%ptr,STAT=err)
                    IF(err/=0) CALL FlagError("Could not allocate equations to solver matrix maps.",err,error,*999)
                    CALL SolverMapping_EquationsToSolverMapsInitialise(solverMapping%equationsSetToSolverMap( &
                      & equationsSetIdx)%equationsToSolverMatrixMapsSm(solverMatrixIdx)% &
                      & linearEquationsToSolverMatrixMaps(equationsMatrixIdx)%ptr,err,error,*999)
                    solverMapping%equationsSetToSolverMap(equationsSetIdx)%equationsToSolverMatrixMapsSm( &
                      & solverMatrixIdx)%linearEquationsToSolverMatrixMaps(equationsMatrixIdx)%ptr% &
                      & equationsMatrixType=SOLVER_MAPPING_EQUATIONS_LINEAR_MATRIX
                    solverMapping%equationsSetToSolverMap(equationsSetIdx)%equationsToSolverMatrixMapsEm( &
                      & equationsMatrixIdx)%equationsMatrixNumber=equationsMatrixIdx
                    solverMapping%equationsSetToSolverMap(equationsSetIdx)%equationsToSolverMatrixMapsEm( &
                      & equationsMatrixIdx)%numberOfSolverMatrices=solverMapping%equationsSetToSolverMap( &
                      & equationsSetIdx)%equationsToSolverMatrixMapsEm(equationsMatrixIdx)%numberOfSolverMatrices+1
                  ENDDO !equationsMatrixIdx
                ELSE
                  ALLOCATE(solverMapping%equationsSetToSolverMap(equationsSetIdx)%equationsToSolverMatrixMapsSm( &
                    & solverMatrixIdx)%linearEquationsToSolverMatrixMaps(0),STAT=err)
                  IF(err/=0) &
                    & CALL FlagError("Could not allocate equations to solver matrix maps sm equations to solver matrix maps.", &
                    & err,error,*999)
                  !Set up linear arrays
                  solverMapping%equationsSetToSolverMap(equationsSetIdx)%equationsToSolverMatrixMapsSm( &
                    & solverMatrixIdx)%numberOfLinearEquationsMatrices=0
                ENDIF
                !Set up nonlinear arrays
                IF(ASSOCIATED(nonlinearMapping)) THEN
                  solverMapping%equationsSetToSolverMap(equationsSetIdx)%equationsToSolverMatrixMapsSm( &
                    & solverMatrixIdx)%numberOfEquationsJacobians=nonlinearMapping%NUMBER_OF_RESIDUAL_VARIABLES
                  ALLOCATE(solverMapping%equationsSetToSolverMap(equationsSetIdx)%equationsToSolverMatrixMapsSm( &
                    & solverMatrixIdx)%jacobianToSolverMatrixMaps(nonlinearMapping%NUMBER_OF_RESIDUAL_VARIABLES),STAT=err)
                  IF(err/=0) CALL FlagError("Could not allocate Jacobian to solver matrix maps.",err,error,*999)
                  DO equationsMatrixIdx=1,nonlinearMapping%NUMBER_OF_RESIDUAL_VARIABLES
                    ALLOCATE(solverMapping%equationsSetToSolverMap(equationsSetIdx)%equationsToSolverMatrixMapsSm( &
                      & solverMatrixIdx)%jacobianToSolverMatrixMaps(equationsMatrixIdx)%ptr,STAT=err)
                    IF(err/=0) CALL FlagError("Could not allocate Jacobian to solver matrix map.",err,error,*999)
                    CALL SolverMapping_JacobianToSolverMapInitialise(solverMapping%equationsSetToSolverMap( &
                      & equationsSetIdx)%equationsToSolverMatrixMapsSm(solverMatrixIdx)% &
                      & jacobianToSolverMatrixMaps(equationsMatrixIdx)%ptr,err,error,*999)
                    solverMapping%equationsSetToSolverMap(equationsSetIdx)%equationsToSolverMatrixMapsJm &
                      & (equationsMatrixIdx)%ptr=>solverMapping%equationsSetToSolverMap(equationsSetIdx)% &
                      & equationsToSolverMatrixMapsSm(solverMatrixIdx)%jacobianToSolverMatrixMaps( &
                      & equationsMatrixIdx)%ptr
                  ENDDO !equationsMatrixIdx
                ENDIF
              ENDIF
              DO variableIdx=1,numberOfVariables
                IF(ASSOCIATED(dynamicMapping)) THEN
                  variableType=solverMapping%createValuesCache%dynamicVariableType(equationsSetIdx)
                  numberOfDynamicEquationsMatrices=dynamicMapping%VAR_TO_EQUATIONS_MATRICES_MAPS(variableType)% &
                    & NUMBER_OF_EQUATIONS_MATRICES
                ELSE
                  IF(ASSOCIATED(nonlinearMapping)) THEN
                    variableType=solverMapping%createValuesCache%residualVariableTypes(variableIdx,equationsSetIdx)
                  ELSE
                    variableType=solverMapping%createValuesCache%matrixVariableTypes(variableIdx,equationsSetIdx, &
                      & solverMatrixIdx)
                    numberOfLinearEquationsMatrices=linearMapping%VAR_TO_EQUATIONS_MATRICES_MAPS(variableType)% &
                      & NUMBER_OF_EQUATIONS_MATRICES
                  ENDIF
                ENDIF

                dependentVariable=>dependentField%VARIABLE_TYPE_MAP(variableType)%ptr
                colDofsMapping=>dependentVariable%DOMAIN_MAPPING
                CALL BOUNDARY_CONDITIONS_VARIABLE_GET(boundaryConditions,dependentVariable,boundaryConditionsVariable, &
                  & err,error,*999)
                IF(ASSOCIATED(dynamicMapping)) THEN
                  !Allocate dynamic equations to solver matrix maps equations column to solver columns maps
                  DO equationsMatrixIdx=1,numberOfDynamicEquationsMatrices
                    matrixNumber=dynamicMapping%VAR_TO_EQUATIONS_MATRICES_MAPS(variableType)% &
                      & EQUATIONS_MATRIX_NUMBERS(equationsMatrixIdx)
                    solverMapping%equationsSetToSolverMap(equationsSetIdx)%equationsToSolverMatrixMapsSm( &
                      & solverMatrixIdx)%dynamicEquationsToSolverMatrixMaps(equationsMatrixIdx)%ptr% &
                      & solverMatrixNumber=solverMatrixIdx
                    solverMapping%equationsSetToSolverMap(equationsSetIdx)%equationsToSolverMatrixMapsSm( &
                      & solverMatrixIdx)%dynamicEquationsToSolverMatrixMaps(equationsMatrixIdx)%ptr% &
                      & equationsMatrixNumber=matrixNumber
                    solverMapping%equationsSetToSolverMap(equationsSetIdx)%equationsToSolverMatrixMapsSm( &
                      & solverMatrixIdx)%dynamicEquationsToSolverMatrixMaps(equationsMatrixIdx)%ptr% &
                      & equationsMatrix=>dynamicMapping%EQUATIONS_MATRIX_TO_VAR_MAPS(matrixNumber)%EQUATIONS_MATRIX
                    numberOfEquationsColumns=dynamicMapping%EQUATIONS_MATRIX_TO_VAR_MAPS(matrixNumber)%NUMBER_OF_COLUMNS
                    ALLOCATE(solverMapping%equationsSetToSolverMap(equationsSetIdx)% &
                      & equationsToSolverMatrixMapsSm(solverMatrixIdx)%dynamicEquationsToSolverMatrixMaps( &
                      & equationsMatrixIdx)%ptr%equationsColToSolverColsMap(numberOfEquationsColumns),STAT=err)
                    IF(err/=0) CALL FlagError("Could not allocate dynamic equations column to solver columns map.", &
                      & err,error,*999)
                  ENDDO !equationsMatrixIdx
                ELSE
                  !Allocate linear equations to solver matrix maps equations column to solver columns maps
                  IF(ASSOCIATED(linearMapping)) THEN
                    DO equationsMatrixIdx=1,numberOfLinearEquationsMatrices
                      matrixNumber=linearMapping%VAR_TO_EQUATIONS_MATRICES_MAPS(variableType)% &
                        & EQUATIONS_MATRIX_NUMBERS(equationsMatrixIdx)
                      solverMapping%equationsSetToSolverMap(equationsSetIdx)%equationsToSolverMatrixMapsSm( &
                        & solverMatrixIdx)%linearEquationsToSolverMatrixMaps(equationsMatrixIdx)%ptr% &
                        & solverMatrixNumber=solverMatrixIdx
                      solverMapping%equationsSetToSolverMap(equationsSetIdx)%equationsToSolverMatrixMapsSm( &
                        & solverMatrixIdx)%linearEquationsToSolverMatrixMaps(equationsMatrixIdx)%ptr% &
                        & equationsMatrixNumber=matrixNumber
                      solverMapping%equationsSetToSolverMap(equationsSetIdx)%equationsToSolverMatrixMapsSm( &
                        & solverMatrixIdx)%linearEquationsToSolverMatrixMaps(equationsMatrixIdx)%ptr% &
                        & equationsMatrix=>linearMapping%EQUATIONS_MATRIX_TO_VAR_MAPS(matrixNumber)%EQUATIONS_MATRIX
                      numberOfEquationsColumns=linearMapping%EQUATIONS_MATRIX_TO_VAR_MAPS(matrixNumber)%NUMBER_OF_COLUMNS
                      ALLOCATE(solverMapping%equationsSetToSolverMap(equationsSetIdx)% &
                        & equationsToSolverMatrixMapsSm(solverMatrixIdx)%linearEquationsToSolverMatrixMaps( &
                        & equationsMatrixIdx)%ptr%equationsColToSolverColsMap(numberOfEquationsColumns),STAT=err)
                      IF(err/=0) CALL FlagError("Could not allocate linear equations column to solver columns map.", &
                        & err,error,*999)
                    ENDDO !equationsMatrixIdx
                  ENDIF
                ENDIF
              ENDDO !variableIdx
              IF(ASSOCIATED(nonlinearMapping)) THEN
                DO equationsMatrixIdx=1,nonlinearMapping%NUMBER_OF_RESIDUAL_VARIABLES
                  solverMapping%equationsSetToSolverMap(equationsSetIdx)%equationsToSolverMatrixMapsSm( &
                    & solverMatrixIdx)%jacobianToSolverMatrixMaps(equationsMatrixIdx)%ptr%solverMatrixNumber= &
                    & solverMatrixIdx
                  solverMapping%equationsSetToSolverMap(equationsSetIdx)%equationsToSolverMatrixMapsSm( &
                    & solverMatrixIdx)%jacobianToSolverMatrixMaps(equationsMatrixIdx)%ptr%jacobianMatrix=> &
                    & nonlinearMapping%JACOBIAN_TO_VAR_MAP(equationsMatrixIdx)%JACOBIAN
                  numberOfEquationsColumns=nonlinearMapping%JACOBIAN_TO_VAR_MAP(equationsMatrixIdx)%NUMBER_OF_COLUMNS
                  ALLOCATE(solverMapping%equationsSetToSolverMap(equationsSetIdx)%equationsToSolverMatrixMapsSm( &
                    & solverMatrixIdx)%jacobianToSolverMatrixMaps(equationsMatrixIdx)%ptr% &
                    & jacobianColToSolverColsMap(numberOfEquationsColumns),STAT=err)
                  IF(err/=0) CALL FlagError("Could not allocate Jacobian column to solver columns map.",err,error,*999)
                ENDDO
              ENDIF
            ENDDO !equationsSetIdx
            
            DO interfaceConditionIdx=1,solverMapping%numberOfInterfaceConditions
              
              interfaceCondition=>solverMapping%interfaceConditions(interfaceConditionIdx)%ptr
              interfaceEquations=>interfaceCondition%INTERFACE_EQUATIONS
              interfaceMapping=>interfaceEquations%INTERFACE_MAPPING

              SELECT CASE(interfaceCondition%method)
              CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD,INTERFACE_CONDITION_PENALTY_METHOD)
                
                !Initialise solver columns to interface condition map
                CALL SolverMapping_SolverColToInterfaceMapInitialise(solverMapping%solverColToEquationsColsMap( &
                  & solverMatrixIdx)%solverColToInterfaceMaps(interfaceConditionIdx),err,error,*999)
                
                !Allocate the solver columns to equations set map arrays
                ALLOCATE(solverMapping%solverColToEquationsColsMap(solverMatrixIdx)%solverColToInterfaceMaps( &
                  & interfaceConditionIdx)%solverColToInterfaceEquationsMaps(numberOfColumns),STAT=err)
                IF(err/=0) CALL FlagError("Could not allocate solver columns to interface equations map.",err,error,*999)
                                
                !Allocate the interface to solver matrix maps sm interface to solver maps
                ALLOCATE(solverMapping%interfaceConditionToSolverMap(interfaceConditionIdx)% &
                  & interfaceToSolverMatrixMapsSm(solverMatrixIdx)%interfaceEquationsToSolverMatrixMaps( &
                  & interfaceMapping%NUMBER_OF_INTERFACE_MATRICES),STAT=err)
                IF(err/=0) CALL FlagError("Could not allocate interface to solver matrix maps sm interface equations "// &
                  & "to solver matrix maps.",err,error,*999)
                
                !Set up interface arrays
                solverMapping%interfaceConditionToSolverMap(interfaceConditionIdx)%interfaceToSolverMatrixMapsSm( &
                  & solverMatrixIdx)%numberOfInterfaceMatrices=interfaceMapping%NUMBER_OF_INTERFACE_MATRICES                
                DO interfaceMatrixIdx=1,interfaceMapping%NUMBER_OF_INTERFACE_MATRICES
                  NULLIFY(solverMapping%interfaceConditionToSolverMap(interfaceConditionIdx)% &
                    & interfaceToSolverMatrixMapsSm(solverMatrixIdx)%interfaceEquationsToSolverMatrixMaps( &
                    & interfaceMatrixIdx)%ptr)
                  ALLOCATE(solverMapping%interfaceConditionToSolverMap(interfaceConditionIdx)% &
                    & interfaceToSolverMatrixMapsSm(solverMatrixIdx)%interfaceEquationsToSolverMatrixMaps( &
                    & interfaceMatrixIdx)%ptr,STAT=err)
                  IF(err/=0) CALL FlagError("Could not allocate interface to solver matrix maps.",err,error,*999)
                  CALL SolverMapping_InterfaceConditionToSolverMapsInitialise(solverMapping% &
                    & interfaceConditionToSolverMap(interfaceConditionIdx)%interfaceToSolverMatrixMapsSm( &
                    & solverMatrixIdx)%interfaceEquationsToSolverMatrixMaps(interfaceMatrixIdx)%ptr,err,error,*999)
                  solverMapping%interfaceConditionToSolverMap(interfaceConditionIdx)%interfaceToSolverMatrixMapsIm( &
                    & interfaceMatrixIdx)%interfaceMatrixNumber=interfaceMatrixIdx
                  solverMapping%interfaceConditionToSolverMap(interfaceConditionIdx)%interfaceToSolverMatrixMapsIm( &
                    & interfaceMatrixIdx)%numberOfSolverMatrices=1
                  
                  dependentVariable=>interfaceMapping%INTERFACE_MATRIX_ROWS_TO_VAR_MAPS(interfaceMatrixIdx)%variable
                  solverMapping%interfaceConditionToSolverMap(interfaceConditionIdx)%interfaceToSolverMatrixMapsSm( &
                    & solverMatrixIdx)%dependentVariables(interfaceMatrixIdx)%ptr=>dependentVariable
                  solverMapping%interfaceConditionToSolverMap(interfaceConditionIdx)%interfaceToSolverMatrixMapsSm( &
                    & solverMatrixIdx)%interfaceEquationsToSolverMatrixMaps(interfaceMatrixIdx)%ptr% &
                    & solverMatrixNumber=solverMatrixIdx
                  solverMapping%interfaceConditionToSolverMap(interfaceConditionIdx)%interfaceToSolverMatrixMapsSm( &
                    & solverMatrixIdx)%interfaceEquationsToSolverMatrixMaps(interfaceMatrixIdx)%ptr% &
                    & interfaceMatrixNumber=interfaceMatrixIdx
                  solverMapping%interfaceConditionToSolverMap(interfaceConditionIdx)%interfaceToSolverMatrixMapsSm( &
                    & solverMatrixIdx)%interfaceEquationsToSolverMatrixMaps(interfaceMatrixIdx)%ptr% &
                    & interfaceMatrix=>interfaceMapping%INTERFACE_MATRIX_ROWS_TO_VAR_MAPS(interfaceMatrixIdx)%INTERFACE_MATRIX
                  numberOfInterfaceRows=interfaceMapping%INTERFACE_MATRIX_ROWS_TO_VAR_MAPS(interfaceMatrixIdx)%NUMBER_OF_ROWS
                  ALLOCATE(solverMapping%interfaceConditionToSolverMap(interfaceConditionIdx)% &
                    & interfaceToSolverMatrixMapsSm(solverMatrixIdx)%interfaceEquationsToSolverMatrixMaps( &
                    & interfaceMatrixIdx)%ptr%interfaceRowToSolverColsMap(numberOfInterfaceRows),STAT=err)
                  IF(err/=0) CALL FlagError("Could not allocate interface column to solver columns map.",err,error,*999)
                ENDDO !interfaceMatrixIdx
                numberOfInterfaceColumns=interfaceMapping%NUMBER_OF_COLUMNS
                ALLOCATE(solverMapping%interfaceConditionToSolverMap(interfaceConditionIdx)% &
                  & interfaceToSolverMatrixMapsSm(solverMatrixIdx)%interfaceColToSolverColsMap( &
                  & numberOfInterfaceColumns),STAT=err)
                IF(err/=0) CALL FlagError("Could not allocate interface column to solver columns map.",err,error,*999)
              CASE(INTERFACE_CONDITION_AUGMENTED_LAGRANGE_METHOD)
                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
              CASE(INTERFACE_CONDITION_POINT_TO_POINT_METHOD)
                CALL FlagError("Not implemented.",err,error,*999)
              CASE DEFAULT
                localError="The interface condition method of "//TRIM(NumberToVString(interfaceCondition%method,"*",err,error))// &
                  & " is invalid."
                CALL FlagError(localError,err,error,*999)
              END SELECT
            ENDDO !interfaceConditionIdx
            
            !Loop over the ranks to ensure that the lowest ranks have the lowest numbered solver variables

            !Allocate dof map to record column reordering
            ALLOCATE(dofMap(solverMapping%columnVariablesList(solverMatrixIdx)%numberOfVariables),STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate dof map.",err,error,*999)
            DO solverVariableIdx=1,solverMapping%columnVariablesList(solverMatrixIdx)%numberOfVariables
              ALLOCATE(dofMap(solverVariableIdx)%ptr(solverMapping%columnVariablesList(solverMatrixIdx)% &
                & variables(solverVariableIdx)%variable%NUMBER_OF_GLOBAL_DOFS),STAT=err)
              IF(err/=0) CALL FlagError("Could not allocate dof map global dof map.",err,error,*999)
              dofMap(solverVariableIdx)%ptr=0
            ENDDO !solverVariableIdx

            ALLOCATE(solverLocalDof(0:COMPUTATIONAL_ENVIRONMENT%NUMBER_COMPUTATIONAL_NODES-1),STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate solver local dof array.",err,error,*999)

            !
            ! 4d Now calculate the solver mappings for each column in rank order
            !
            
            numberOfGlobalSolverDofs=0
            solverGlobalDof=0
            solverLocalDof=0
            DO dofType=1,2
              DO rank=0,COMPUTATIONAL_ENVIRONMENT%NUMBER_COMPUTATIONAL_NODES-1
                
                DO solverVariableIdx=1,solverMapping%columnVariablesList(solverMatrixIdx)%numberOfVariables
                  
                  IF (solverMapping%numberOfInterfaceConditions>0) THEN
                    ! Ensure that the dof_offset is calculated as a sum of the number of dofs in the diagonal entries of the solver
                    ! matrix (ie the sum of the number of solver dofs in each equation set).
                    ! Note that this may not work when running problems in parallel, however, note that interfaces can not currently
                    ! be used in parallel either, and the following code only executes if there are interface conditions present.
                    tempOffset = 0
                    DO solverVariableIdxTemp=1,solverVariableIdx
                      DO globalDof=1,SIZE(dofMap(solverVariableIdxTemp)%ptr)
                        IF (dofMap(solverVariableIdxTemp)%PTR(globalDof)>0) THEN
                          tempOffset=tempOffset+1
                        ENDIF
                      ENDDO
                    ENDDO
                    globalDofsOffset=tempOffset
                    localDofsOffset=tempOffset
                  ELSE
                    globalDofsOffset=solverGlobalDof
                    localDofsOffset=solverLocalDof(rank)
                  ENDIF
                  
                  variableType=solverMapping%columnVariablesList(solverMatrixIdx)%variables(solverVariableIdx)%variableType
                  
                  DO equationsIdx=1,solverMapping%numberOfEquationsSets+solverMapping%numberOfInterfaceConditions
                    
                    !Get columns list
                    CALL List_Sort(rankGlobalColsLists(dofType,equationsIdx,solverVariableIdx,rank)%ptr,err,error,*999)
                    CALL List_DetachAndDestroy(rankGlobalColsLists(dofType,equationsIdx,solverVariableIdx,rank)%ptr, &
                      & numberOfRankCols,rankGlobalColsList,err,error,*999)
                    
                    IF(numberOfRankCols>0) THEN

                      solverGlobalDof=globalDofsOffset
                      solverLocalDof(rank)=localDofsOffset
                    
                      equationType=subMatrixInformation(1,equationsIdx,solverVariableIdx)
                      SELECT CASE(equationType)
                      CASE(SOLVER_MAPPING_EQUATIONS_EQUATIONS_SET)
                        
                        equationsSetIdx=subMatrixInformation(2,equationsIdx,solverVariableIdx)
                       
                        !The pointers below have been checked for association above.
                        equationsSet=>solverMapping%equationsSets(equationsSetIdx)%ptr
                        equations=>equationsSet%equations
                        equationsMapping=>equations%EQUATIONS_MAPPING
                        dynamicMapping=>equationsMapping%DYNAMIC_MAPPING
                        linearMapping=>equationsMapping%LINEAR_MAPPING
                        nonlinearMapping=>equationsMapping%NONLINEAR_MAPPING
                        dependentField=>equationsSet%dependent%DEPENDENT_FIELD
 
                        numberOfDynamicEquationsMatrices=0
                        numberOfLinearEquationsMatrices=0
                        IF(ASSOCIATED(dynamicMapping)) numberOfDynamicEquationsMatrices=dynamicMapping% &
                          & VAR_TO_EQUATIONS_MATRICES_MAPS(variableType)%NUMBER_OF_EQUATIONS_MATRICES
                        IF(ASSOCIATED(linearMapping)) numberOfLinearEquationsMatrices=linearMapping% &
                          & VAR_TO_EQUATIONS_MATRICES_MAPS(variableType)%NUMBER_OF_EQUATIONS_MATRICES
                        
                        !Loop over the variables
                        
                        dependentVariable=>solverMapping%columnVariablesList(solverMatrixIdx)%variables(solverVariableIdx)% &
                          & variable
                        colDofsMapping=>dependentVariable%DOMAIN_MAPPING
                        CALL BOUNDARY_CONDITIONS_VARIABLE_GET(boundaryConditions,dependentVariable,boundaryConditionsVariable, &
                          & err,error,*999)
                        
                        DO globalDofIdx=1,numberOfRankCols
                          globalDof=rankGlobalColsList(1,globalDofIdx)
                          localDof=rankGlobalColsList(2,globalDofIdx)
                          !dofType=rankGlobalColsList(3,globalDofIdx)
                          includeColumn=rankGlobalColsList(3,globalDofIdx)==1                      
                          constrainedDof=rankGlobalColsList(3,globalDofIdx)==2                      
                          variableIdx=rankGlobalColsList(4,globalDofIdx)
                          globalDofCouplingNumber=rankGlobalColsList(5,globalDofIdx)
                          IF(globalDofCouplingNumber>0) THEN
                            colEquationCols=>columnCouplings%dofCouplings(globalDofCouplingNumber)%ptr
                            IF(ASSOCIATED(colEquationCols)) THEN
                              numberColEquationsCols=colEquationCols%numberOfDofs
                            ELSE
                              CALL FlagError("Dof coupling is not associated for global dof coupling number "// &
                                & TRIM(NumberToVstring(globalDofCouplingNumber,"*",err,error))//".",err,error,*999)
                            END IF
                          ELSE
                            numberColEquationsCols=1
                            dummyDofCoupling%globalDofs(1)=globalDof
                            dummyDofCoupling%localDofs(1)=localDof
                            dummyDofCoupling%coefficients(1)=1.0_DP
                            colEquationCols=>dummyDofCoupling
                          END IF
                        
                          IF(includeColumn) THEN
                            !DOF is not fixed so map the variable/equation dof to a new solver dof
                            
                            IF(dofType==2) THEN
                              solverGlobalDof=dofMap(solverVariableIdx)%ptr(globalDof)
                            ELSE
                              solverGlobalDof=solverGlobalDof+1
                              dofMap(solverVariableIdx)%ptr(globalDof)=solverGlobalDof
                            ENDIF
                            
                            solverLocalDof(rank)=solverLocalDof(rank)+1
                            
                            IF(rank==myrank) THEN
                              
                              IF(.NOT.variableRankProcessed(solverVariableIdx,rank)) THEN
                                
                                !Set up the column domain mappings.
                                CALL DOMAIN_MAPPINGS_MAPPING_GLOBAL_INITIALISE(colDomainMapping%GLOBAL_TO_LOCAL_MAP( &
                                  & solverGlobalDof),err,error,*999)
                                !There are no ghosted cols for the solver matrices so there is only 1 domain for the global to
                                !local map.
                                !Allocate the global to local map arrays
                                ALLOCATE(colDomainMapping%GLOBAL_TO_LOCAL_MAP(solverGlobalDof)%LOCAL_NUMBER(1),STAT=err)
                                IF(err/=0) CALL FlagError("Could not allocate column domain global to local map local number.", &
                                  & err,error,*999)
                                ALLOCATE(colDomainMapping%GLOBAL_TO_LOCAL_MAP(solverGlobalDof)%DOMAIN_NUMBER(1),STAT=err)
                                IF(err/=0) CALL FlagError("Could not allocate column domain global to local map domain number.", &
                                  & err,error,*999)
                                ALLOCATE(colDomainMapping%GLOBAL_TO_LOCAL_MAP(solverGlobalDof)%LOCAL_TYPE(1),STAT=err)
                                IF(err/=0) CALL FlagError("Could not allocate column domain global to local map domain number.", &
                                  & err,error,*999)
                                !Set up the global to local mappings
                                colDomainMapping%GLOBAL_TO_LOCAL_MAP(solverGlobalDof)%NUMBER_OF_DOMAINS=1
                                colDomainMapping%GLOBAL_TO_LOCAL_MAP(solverGlobalDof)%LOCAL_NUMBER(1)=solverLocalDof(rank)
                                colDomainMapping%GLOBAL_TO_LOCAL_MAP(solverGlobalDof)%DOMAIN_NUMBER(1)=rank
                                colDomainMapping%GLOBAL_TO_LOCAL_MAP(solverGlobalDof)%LOCAL_TYPE(1)=DOMAIN_LOCAL_INTERNAL
                                
                                !Set up the solver column -> equations column mappings. 1-1 as no coupling yet
!!TODO
                                !Set up the solver dofs -> variable dofs map
                                !Initialise
                                CALL SolverMapping_SolverDofToVariableMapInitialise(solverMapping% &
                                  & solverColToEquationsColsMap(solverMatrixIdx)%solverDofToVariableMaps( &
                                  & solverLocalDof(rank)),err,error,*999)
                                !Allocate the solver dofs to variable dofs arrays
!!TODO: allow for multiple equations sets in the column
                                !No coupling so there is only one equations set at the moment
                                ALLOCATE(solverMapping%solverColToEquationsColsMap(solverMatrixIdx)% &
                                  & solverDofToVariableMaps(solverLocalDof(rank))%equationsTypes(numberColEquationsCols),STAT=err)
                                IF(err/=0) CALL FlagError("Could not allocate solver dof to variable maps equations types.", &
                                  & err,error,*999)
                                ALLOCATE(solverMapping%solverColToEquationsColsMap(solverMatrixIdx)% &
                                  & solverDofToVariableMaps(solverLocalDof(rank))%equationsIndices(numberColEquationsCols),STAT=err)
                                IF(err/=0) CALL FlagError("Could not allocate solver dof to variable maps equations indices.", &
                                  & err,error,*999)
                                ALLOCATE(solverMapping%solverColToEquationsColsMap(solverMatrixIdx)% &
                                  & solverDofToVariableMaps(solverLocalDof(rank))%variable(numberColEquationsCols),STAT=err)
                                IF(err/=0) CALL FlagError("Could not allocate solver dof to variable maps variable type.", &
                                  & err,error,*999)
                                ALLOCATE(solverMapping%solverColToEquationsColsMap(solverMatrixIdx)% &
                                  & solverDofToVariableMaps(solverLocalDof(rank))%variableDof(numberColEquationsCols),STAT=err)
                                IF(err/=0) CALL FlagError("Could not allocate solver dof to variable maps variable dof.", &
                                  & err,error,*999)
                                ALLOCATE(solverMapping%solverColToEquationsColsMap(solverMatrixIdx)% &
                                  & solverDofToVariableMaps(solverLocalDof(rank))%variableCoefficient(numberColEquationsCols), &
                                  & STAT=err)
                                IF(err/=0) CALL FlagError("Could not allocate solver dof to variable maps variable coefficient.", &
                                  & err,error,*999)
                                ALLOCATE(solverMapping%solverColToEquationsColsMap(solverMatrixIdx)% &
                                  & solverDofToVariableMaps(solverLocalDof(rank))%additiveConstant(numberColEquationsCols),STAT=err)
                                IF(err/=0) CALL FlagError("Could not allocate solver dof to variable maps additive constant.", &
                                  & err,error,*999)
                                !Set the solver -> equations mappings
                                solverMapping%solverColToEquationsColsMap(solverMatrixIdx)%solverDofToVariableMaps( &
                                  & solverLocalDof(rank))%numberOfEquations=numberColEquationsCols
                                DO colEquationsColIdx=1,numberColEquationsCols
                                  solverMapping%solverColToEquationsColsMap(solverMatrixIdx)%solverDofToVariableMaps( &
                                    & solverLocalDof(rank))%equationsTypes(colEquationsColIdx)= &
                                    & SOLVER_MAPPING_EQUATIONS_EQUATIONS_SET
                                  solverMapping%solverColToEquationsColsMap(solverMatrixIdx)%solverDofToVariableMaps( &
                                    & solverLocalDof(rank))%equationsIndices(colEquationsColIdx)=equationsSetIdx
                                  solverMapping%solverColToEquationsColsMap(solverMatrixIdx)%solverDofToVariableMaps( &
                                    & solverLocalDof(rank))%variable(colEquationsColIdx)%ptr=>dependentVariable
                                  solverMapping%solverColToEquationsColsMap(solverMatrixIdx)%solverDofToVariableMaps( &
                                    & solverLocalDof(rank))%variableDof(colEquationsColIdx)= &
                                    & colEquationCols%localDofs(colEquationsColIdx)
                                  solverMapping%solverColToEquationsColsMap(solverMatrixIdx)%solverDofToVariableMaps( &
                                    & solverLocalDof(rank))%variableCoefficient(colEquationsColIdx)= &
                                    & colEquationCols%coefficients(colEquationsColIdx)
                                  solverMapping%solverColToEquationsColsMap(solverMatrixIdx)%solverDofToVariableMaps( &
                                    & solverLocalDof(rank))%additiveConstant(colEquationsColIdx)=0.0_DP
                                ENDDO !colEquationsColIdx
                              ENDIF
                              !Set up the equations variables -> solver columns mapping
                              DO colEquationsColIdx=1,numberColEquationsCols
                                eqnLocalDof=colEquationCols%localDofs(colEquationsColIdx)
                                couplingCoefficient=colEquationCols%coefficients(colEquationsColIdx)
                                solverMapping%equationsSetToSolverMap(equationsSetIdx)%equationsToSolverMatrixMapsSm( &
                                  & solverMatrixIdx)%variableToSolverColMaps(variableIdx)%columnNumbers(eqnLocalDof)= &
                                  & solverGlobalDof
                                solverMapping%equationsSetToSolverMap(equationsSetIdx)%equationsToSolverMatrixMapsSm( &
                                  & solverMatrixIdx)%variableToSolverColMaps(variableIdx)%couplingCoefficients( &
                                  & eqnLocalDof)=couplingCoefficient
                                solverMapping%equationsSetToSolverMap(equationsSetIdx)%equationsToSolverMatrixMapsSm( &
                                  & solverMatrixIdx)%variableToSolverColMaps(variableIdx)%additiveConstants( &
                                  & eqnLocalDof)=0.0_DP
                                !Set up the equations columns -> solver columns mapping
                                DO matrixTypeIdx=1,subMatrixList(0,equationsIdx,solverVariableIdx)
                                  matrixType=subMatrixList(matrixTypeIdx,equationsIdx,solverVariableIdx)
                                  SELECT CASE(matrixType)
                                  CASE(SOLVER_MAPPING_EQUATIONS_DYNAMIC_MATRIX)
                                    !Dynamic matrix
                                    DO equationsMatrixIdx=1,numberOfDynamicEquationsMatrices
                                      matrixNumber=dynamicMapping%VAR_TO_EQUATIONS_MATRICES_MAPS(variableType)% &
                                        & EQUATIONS_MATRIX_NUMBERS(equationsMatrixIdx)
                                      equationsColumn=dynamicMapping%VAR_TO_EQUATIONS_MATRICES_MAPS(variableType)% &
                                        & DOF_TO_COLUMNS_MAPS(equationsMatrixIdx)%COLUMN_DOF(localDof)
                                      !Allocate the equation to solver map column items.
                                      ALLOCATE(solverMapping%equationsSetToSolverMap(equationsSetIdx)% &
                                        & equationsToSolverMatrixMapsSm(solverMatrixIdx)% &
                                        & dynamicEquationsToSolverMatrixMaps(equationsMatrixIdx)%ptr% &
                                        & equationsColToSolverColsMap(equationsColumn)%solverCols(1), &
                                      & STAT=err)
                                      IF(err/=0) &
                                        & CALL  FlagError("Could not allocate dynamic equations column to solver columns map "// &
                                        & "solver columns.",err,error,*999)
                                      ALLOCATE(solverMapping%equationsSetToSolverMap(equationsSetIdx)% &
                                        & equationsToSolverMatrixMapsSm(solverMatrixIdx)% &
                                        & dynamicEquationsToSolverMatrixMaps(equationsMatrixIdx)%ptr% &
                                        & equationsColToSolverColsMap(equationsColumn)%couplingCoefficients(1),STAT=err)
                                      IF(err/=0) &
                                        & CALL FlagError("Could not allocate dynamic equations column to solver columns map "// &
                                        & "coupling coefficients.",err,error,*999)
                                      solverMapping%equationsSetToSolverMap(equationsSetIdx)% &
                                        & equationsToSolverMatrixMapsSm(solverMatrixIdx)% &
                                        & dynamicEquationsToSolverMatrixMaps(equationsMatrixIdx)%ptr% &
                                        & equationsColToSolverColsMap(equationsColumn)%numberOfSolverCols=1
                                      solverMapping%equationsSetToSolverMap(equationsSetIdx)% &
                                        & equationsToSolverMatrixMapsSm(solverMatrixIdx)% &
                                        & dynamicEquationsToSolverMatrixMaps(equationsMatrixIdx)%ptr% &
                                        & equationsColToSolverColsMap(equationsColumn)%solverCols(1)=solverGlobalDof
                                      solverMapping%equationsSetToSolverMap(equationsSetIdx)% &
                                        & equationsToSolverMatrixMapsSm(solverMatrixIdx)% &
                                        & dynamicEquationsToSolverMatrixMaps(equationsMatrixIdx)%ptr% &
                                        & equationsColToSolverColsMap(equationsColumn)%couplingCoefficients(1)= &
                                        & dynamicMapping%EQUATIONS_MATRIX_TO_VAR_MAPS(matrixNumber)%MATRIX_COEFFICIENT
                                    ENDDO !equationsMatrixIdx
                                  CASE(SOLVER_MAPPING_EQUATIONS_LINEAR_MATRIX)
                                    DO equationsMatrixIdx=1,numberOfLinearEquationsMatrices
                                      matrixNumber=linearMapping%VAR_TO_EQUATIONS_MATRICES_MAPS(variableType)% &
                                        & EQUATIONS_MATRIX_NUMBERS(equationsMatrixIdx)
                                      equationsColumn=linearMapping%VAR_TO_EQUATIONS_MATRICES_MAPS(variableType)% &
                                        & DOF_TO_COLUMNS_MAPS(equationsMatrixIdx)%COLUMN_DOF(localDof)
                                      !Allocate the equation to solver map column items.
                                      ALLOCATE(solverMapping%equationsSetToSolverMap(equationsSetIdx)% &
                                        & equationsToSolverMatrixMapsSm(solverMatrixIdx)% &
                                        & linearEquationsToSolverMatrixMaps(equationsMatrixIdx)%ptr% &
                                        & equationsColToSolverColsMap(equationsColumn)%solverCols(1), &
                                        & STAT=err)
                                      IF(err/=0) &
                                        & CALL FlagError("Could not allocate linear equations column to solver columns map "// &
                                        & "solver columns.",err,error,*999)
                                      ALLOCATE(solverMapping%equationsSetToSolverMap(equationsSetIdx)% &
                                        & equationsToSolverMatrixMapsSm(solverMatrixIdx)% &
                                        & linearEquationsToSolverMatrixMaps(equationsMatrixIdx)%ptr% &
                                        & equationsColToSolverColsMap(equationsColumn)%couplingCoefficients(1),STAT=err)
                                      IF(err/=0) &
                                        & CALL FlagError("Could not allocate linear equations column to solver columns map "// &
                                        & "coupling coefficients.",err,error,*999)
                                      solverMapping%equationsSetToSolverMap(equationsSetIdx)% &
                                        & equationsToSolverMatrixMapsSm(solverMatrixIdx)% &
                                        & linearEquationsToSolverMatrixMaps(equationsMatrixIdx)%ptr% &
                                        & equationsColToSolverColsMap(equationsColumn)%numberOfSolverCols=1
                                      solverMapping%equationsSetToSolverMap(equationsSetIdx)% &
                                        & equationsToSolverMatrixMapsSm(solverMatrixIdx)% &
                                        & linearEquationsToSolverMatrixMaps(equationsMatrixIdx)%ptr% &
                                        & equationsColToSolverColsMap(equationsColumn)%solverCols(1)=solverGlobalDof
                                      solverMapping%equationsSetToSolverMap(equationsSetIdx)% &
                                        & equationsToSolverMatrixMapsSm(solverMatrixIdx)% &
                                        & linearEquationsToSolverMatrixMaps(equationsMatrixIdx)%ptr% &
                                        & equationsColToSolverColsMap(equationsColumn)%couplingCoefficients(1)= &
                                        & linearMapping%EQUATIONS_MATRIX_TO_VAR_MAPS(matrixNumber)%MATRIX_COEFFICIENT
                                    ENDDO !equationsMatrixIdx
                                  CASE(SOLVER_MAPPING_EQUATIONS_NONLINEAR_MATRIX)
                                    DO equationsMatrixIdx=1,nonlinearMapping%NUMBER_OF_RESIDUAL_VARIABLES
                                      IF(nonlinearMapping%VAR_TO_JACOBIAN_MAP(equationsMatrixIdx)%VARIABLE_TYPE==variableType) EXIT
                                    ENDDO !equationsMatrixIdx
                                    jacobianColumn=nonlinearMapping%VAR_TO_JACOBIAN_MAP(equationsMatrixIdx)% &
                                      & DOF_TO_COLUMNS_MAP(localDof)
                                    !Allocate the Jacobian to solver map column items.
                                    ALLOCATE(solverMapping%equationsSetToSolverMap(equationsSetIdx)% &
                                      & equationsToSolverMatrixMapsSm(solverMatrixIdx)%jacobianToSolverMatrixMaps( &
                                      & equationsMatrixIdx)%ptr%jacobianColToSolverColsMap(jacobianColumn)%solverCols(1), &
                                      & STAT=err)
                                    IF(err/=0) &
                                      & CALL FlagError("Could not allocate Jacobian column to solver columns map solver columns.", &
                                      & err,error,*999)
                                    ALLOCATE(solverMapping%equationsSetToSolverMap(equationsSetIdx)% &
                                      & equationsToSolverMatrixMapsSm(solverMatrixIdx)%jacobianToSolverMatrixMaps( &
                                      & equationsMatrixIdx)%ptr%jacobianColToSolverColsMap(jacobianColumn)% &
                                      & couplingCoefficients(1),STAT=err)
                                    IF(err/=0) THEN
                                      localError="Could not allocate Jacobain column to solver columns map coupling coefficients."
                                      CALL FlagError(localError,err,error,*999)
                                    ENDIF
                                    solverMapping%equationsSetToSolverMap(equationsSetIdx)% &
                                      & equationsToSolverMatrixMapsSm(solverMatrixIdx)%jacobianToSolverMatrixMaps( &
                                      & equationsMatrixIdx)%ptr%jacobianColToSolverColsMap(jacobianColumn)% &
                                      & numberOfSolverCols=1
                                    solverMapping%equationsSetToSolverMap(equationsSetIdx)% &
                                      & equationsToSolverMatrixMapsSm(solverMatrixIdx)%jacobianToSolverMatrixMaps( &
                                      & equationsMatrixIdx)%ptr%jacobianColToSolverColsMap(jacobianColumn)%solverCols(1)= &
                                      & solverGlobalDof
                                    solverMapping%equationsSetToSolverMap(equationsSetIdx)% &
                                      & equationsToSolverMatrixMapsSm(solverMatrixIdx)%jacobianToSolverMatrixMaps( &
                                      & equationsMatrixIdx)%ptr%jacobianColToSolverColsMap(jacobianColumn)% &
                                      & couplingCoefficients(1)=nonlinearMapping%JACOBIAN_TO_VAR_MAP(equationsMatrixIdx)% &
                                      & JACOBIAN_COEFFICIENT
                                  CASE DEFAULT
                                    localError="The equations matrix type of "// &
                                      & TRIM(NumberToVString(matrixType,"*",err,error))// &
                                      & " is invalid."
                                    CALL FlagError(localError,err,error,*999)
                                  END SELECT
                                ENDDO !matrixTypeIdx
                              END DO !colEquationsColIdx
                            ELSE !rank /= myrank
                              
                              IF(.NOT.variableRankProcessed(solverVariableIdx,rank)) THEN
                                
                                !Set up the column domain mappings.
                                CALL DOMAIN_MAPPINGS_MAPPING_GLOBAL_INITIALISE(colDomainMapping%GLOBAL_TO_LOCAL_MAP( &
                                  & solverGlobalDof),err,error,*999)
                                !There are no ghosted cols for the solver matrices so there is only 1 domain for the global to
                                !local map.
                                !Allocate the global to local map arrays
                                ALLOCATE(colDomainMapping%GLOBAL_TO_LOCAL_MAP(solverGlobalDof)%LOCAL_NUMBER(1),STAT=err)
                                IF(err/=0) CALL FlagError("Could not allocate column domain global to local map local number.", &
                                  & err,error,*999)
                                ALLOCATE(colDomainMapping%GLOBAL_TO_LOCAL_MAP(solverGlobalDof)%DOMAIN_NUMBER(1),STAT=err)
                                IF(err/=0) CALL FlagError("Could not allocate column domain global to local map domain number.", &
                                  & err,error,*999)
                                ALLOCATE(colDomainMapping%GLOBAL_TO_LOCAL_MAP(solverGlobalDof)%LOCAL_TYPE(1),STAT=err)
                                IF(err/=0) CALL FlagError("Could not allocate column domain global to local map domain number.", &
                                  & err,error,*999)
                                !Set up the global to local mappings
                                colDomainMapping%GLOBAL_TO_LOCAL_MAP(solverGlobalDof)%NUMBER_OF_DOMAINS=1
                                colDomainMapping%GLOBAL_TO_LOCAL_MAP(solverGlobalDof)%LOCAL_NUMBER(1)=solverLocalDof(rank)
                                colDomainMapping%GLOBAL_TO_LOCAL_MAP(solverGlobalDof)%DOMAIN_NUMBER(1)=rank
                                colDomainMapping%GLOBAL_TO_LOCAL_MAP(solverGlobalDof)%LOCAL_TYPE(1)=DOMAIN_LOCAL_INTERNAL

                              ENDIF
                              
                            ENDIF !rank == myrank
                          ELSE IF(constrainedDof) THEN
                            !Do nothing, this is set up above
                          ELSE
                            IF(rank==myrank) THEN
                              !Set up the equations variables -> solver columns mapping
                              solverMapping%equationsSetToSolverMap(equationsSetIdx)%equationsToSolverMatrixMapsSm( &
                                & solverMatrixIdx)%variableToSolverColMaps(variableIdx)%columnNumbers(localDof)=0
                              solverMapping%equationsSetToSolverMap(equationsSetIdx)%equationsToSolverMatrixMapsSm( &
                                & solverMatrixIdx)%variableToSolverColMaps(variableIdx)%couplingCoefficients(localDof)=0.0_DP
                              solverMapping%equationsSetToSolverMap(equationsSetIdx)%equationsToSolverMatrixMapsSm( &
                                & solverMatrixIdx)%variableToSolverColMaps(variableIdx)%additiveConstants(localDof)=0.0_DP
                              DO matrixTypeIdx=1,subMatrixList(0,equationsIdx,solverVariableIdx)
                                matrixType=subMatrixList(matrixTypeIdx,equationsIdx,solverVariableIdx)
                                SELECT CASE(matrixType)
                                CASE(SOLVER_MAPPING_EQUATIONS_DYNAMIC_MATRIX)
                                  !Set up the equations columns -> solver columns mapping
                                  DO equationsMatrixIdx=1,numberOfDynamicEquationsMatrices
                                    matrixNumber=dynamicMapping%VAR_TO_EQUATIONS_MATRICES_MAPS(variableType)% &
                                      & EQUATIONS_MATRIX_NUMBERS(equationsMatrixIdx)
                                    equationsColumn=dynamicMapping%VAR_TO_EQUATIONS_MATRICES_MAPS(variableType)% &
                                      & DOF_TO_COLUMNS_MAPS(equationsMatrixIdx)%COLUMN_DOF(localDof)
                                    solverMapping%equationsSetToSolverMap(equationsSetIdx)% &
                                      & equationsToSolverMatrixMapsSm(solverMatrixIdx)% &
                                      & dynamicEquationsToSolverMatrixMaps(equationsMatrixIdx)%ptr% &
                                      & equationsColToSolverColsMap(equationsColumn)%numberOfSolverCols=0
                                  ENDDO !equationsMatrixIdx
                                CASE(SOLVER_MAPPING_EQUATIONS_LINEAR_MATRIX)
                                  !Set up the equations columns -> solver columns mapping
                                  DO equationsMatrixIdx=1,numberOfLinearEquationsMatrices
                                    matrixNumber=linearMapping%VAR_TO_EQUATIONS_MATRICES_MAPS(variableType)% &
                                      & EQUATIONS_MATRIX_NUMBERS(equationsMatrixIdx)
                                    equationsColumn=linearMapping%VAR_TO_EQUATIONS_MATRICES_MAPS(variableType)% &
                                      & DOF_TO_COLUMNS_MAPS(equationsMatrixIdx)%COLUMN_DOF(localDof)
                                    solverMapping%equationsSetToSolverMap(equationsSetIdx)% &
                                      & equationsToSolverMatrixMapsSm(solverMatrixIdx)% &
                                      & linearEquationsToSolverMatrixMaps(equationsMatrixIdx)%ptr% &
                                      & equationsColToSolverColsMap(equationsColumn)%numberOfSolverCols=0
                                  ENDDO !equationsMatrixIdx
                                CASE(SOLVER_MAPPING_EQUATIONS_NONLINEAR_MATRIX)
                                  DO equationsMatrixIdx=1,nonlinearMapping%NUMBER_OF_RESIDUAL_VARIABLES
                                    IF(nonlinearMapping%VAR_TO_JACOBIAN_MAP(equationsMatrixIdx)%VARIABLE_TYPE==variableType) EXIT
                                  ENDDO !equationsMatrixIdx
                                  jacobianColumn=nonlinearMapping%VAR_TO_JACOBIAN_MAP(equationsMatrixIdx)% &
                                    & DOF_TO_COLUMNS_MAP(localDof)
                                  solverMapping%equationsSetToSolverMap(equationsSetIdx)% &
                                    & equationsToSolverMatrixMapsSm(solverMatrixIdx)%jacobianToSolverMatrixMaps( &
                                    & equationsMatrixIdx)%ptr%jacobianColToSolverColsMap(jacobianColumn)% &
                                    & numberOfSolverCols=0
                                CASE DEFAULT
                                  localError="The equations matrix type of "// &
                                    & TRIM(NumberToVString(matrixType,"*",err,error))//" is invalid."
                                  CALL FlagError(localError,err,error,*999)
                                END SELECT
                              ENDDO !matrixTypeIdx
                            ENDIF !rank == myrank
                          ENDIF !include_column
                        ENDDO !globalDof

                      CASE(SOLVER_MAPPING_EQUATIONS_INTERFACE_CONDITION,SOLVER_MAPPING_EQUATIONS_INTERFACE_TRANSPOSE)
                        
                        !Now handle the interface condition rows and columns
                        
                        interfaceConditionIdx=subMatrixInformation(2,equationsIdx,solverVariableIdx)
                        interfaceMatrixIdx=subMatrixInformation(3,equationsIdx,solverVariableIdx)
                        
                        !The pointers below have been checked for association above.
                        interfaceCondition=>solverMapping%interfaceConditions(interfaceConditionIdx)%ptr
                        
                        SELECT CASE(interfaceCondition%METHOD)
                        CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD,INTERFACE_CONDITION_PENALTY_METHOD)
                          interfaceEquations=>interfaceCondition%INTERFACE_EQUATIONS
                          interfaceMapping=>interfaceEquations%INTERFACE_MAPPING

                          !Loop over the variables
                          ! This is not only a Lagrange variable (it could be a equationset variable) - rename for clarity.
                          
                          lagrangeVariable=>solverMapping%columnVariablesList(solverMatrixIdx)%variables(solverVariableIdx)% &
                            & variable
                         
                          DO globalDofIdx=1,numberOfRankCols
                            globalDof=rankGlobalColsList(1,globalDofIdx)
                            localDof=rankGlobalColsList(2,globalDofIdx)
                            !dofType=rankGlobalColsList(3,globalDofIdx)
                            includeColumn=rankGlobalColsList(3,globalDofIdx)==1
                            constrainedDof=rankGlobalColsList(3,globalDofIdx)==2
                            globalDofCouplingNumber=rankGlobalColsList(5,globalDofIdx)

                            IF(globalDofCouplingNumber>0) THEN
                              colEquationCols=>columnCouplings%dofCouplings(globalDofCouplingNumber)%ptr
                              IF(ASSOCIATED(colEquationCols)) THEN
                                numberColEquationsCols=colEquationCols%numberOfDofs
                              ELSE
                                CALL FlagError("Dof coupling is not associated for global dof coupling number "// &
                                  & TRIM(NumberToVstring(globalDofCouplingNumber,"*",err,error))//".",err,error,*999)
                              END IF
                            ELSE
                              numberColEquationsCols=1
                              dummyDofCoupling%globalDofs(1)=globalDof
                              dummyDofCoupling%localDofs(1)=localDof
                              dummyDofCoupling%coefficients(1)=1.0_DP
                              colEquationCols=>dummyDofCoupling
                            END IF
                             
                            IF(includeColumn) THEN
                              !DOF is not fixed so map the variable/equation dof to a new solver dof
                              
                              IF(dofType==2) THEN
                                !Ghosted, reuse global dof
                                solverGlobalDof=dofMap(solverVariableIdx)%ptr(globalDof)
                              ELSE
                                solverGlobalDof=solverGlobalDof+1
                                dofMap(solverVariableIdx)%ptr(globalDof)=solverGlobalDof
                              ENDIF
                              
                              solverLocalDof(rank)=solverLocalDof(rank)+1
                              
                              IF(rank==myrank) THEN
                                                                 
                                IF(.NOT.variableRankProcessed(solverVariableIdx,rank)) THEN
                                  
                                  !Set up the column domain mappings.
                                  CALL DOMAIN_MAPPINGS_MAPPING_GLOBAL_INITIALISE(colDomainMapping%GLOBAL_TO_LOCAL_MAP( &
                                    & solverGlobalDof),err,error,*999)
                                  !There are no ghosted cols for the solver matrices so there is only 1 domain for the global to
                                  !local map.
                                  !Allocate the global to local map arrays
                                  ALLOCATE(colDomainMapping%GLOBAL_TO_LOCAL_MAP(solverGlobalDof)%LOCAL_NUMBER(1),STAT=err)
                                  IF(err/=0) &
                                    & CALL FlagError("Could not allocate column domain global to local map local number.", &
                                    & err,error,*999)
                                  ALLOCATE(colDomainMapping%GLOBAL_TO_LOCAL_MAP(solverGlobalDof)%DOMAIN_NUMBER(1),STAT=err)
                                  IF(err/=0) &
                                    & CALL FlagError("Could not allocate column domain global to local map domain number.", &
                                    & err,error,*999)
                                  ALLOCATE(colDomainMapping%GLOBAL_TO_LOCAL_MAP(solverGlobalDof)%LOCAL_TYPE(1),STAT=err)
                                  IF(err/=0) &
                                    & CALL FlagError("Could not allocate column domain global to local map domain number.", &
                                    & err,error,*999)
                                  !Set up the global to local mappings
                                  colDomainMapping%GLOBAL_TO_LOCAL_MAP(solverGlobalDof)%NUMBER_OF_DOMAINS=1
                                  colDomainMapping%GLOBAL_TO_LOCAL_MAP(solverGlobalDof)%LOCAL_NUMBER(1)=solverLocalDof(rank)
                                  colDomainMapping%GLOBAL_TO_LOCAL_MAP(solverGlobalDof)%DOMAIN_NUMBER(1)=rank
                                  colDomainMapping%GLOBAL_TO_LOCAL_MAP(solverGlobalDof)%LOCAL_TYPE(1)=DOMAIN_LOCAL_INTERNAL
                                  
                                  !Set up the solver column -> equations column mappings.
                                  !Set up the solver dofs -> variable dofs map
                                  !Initialise
                                  CALL SolverMapping_SolverDofToVariableMapInitialise(solverMapping% &
                                    & solverColToEquationsColsMap(solverMatrixIdx)%solverDofToVariableMaps( &
                                    & solverLocalDof(rank)),err,error,*999)
                                  !Allocate the solver dofs to variable dofs arrays
!!TODO: allow for multiple equations sets in the column
                                  ALLOCATE(solverMapping%solverColToEquationsColsMap(solverMatrixIdx)% &
                                    & solverDofToVariableMaps(solverLocalDof(rank))%equationsTypes(numberColEquationsCols),STAT=err)
                                  IF(err/=0) CALL FlagError("Could not allocate solver dof to variable maps equations types.", &
                                    & err,error,*999)
                                  ALLOCATE(solverMapping%solverColToEquationsColsMap(solverMatrixIdx)% &
                                    & solverDofToVariableMaps(solverLocalDof(rank))%equationsIndices(numberColEquationsCols), &
                                    & STAT=err)
                                  IF(err/=0) CALL FlagError("Could not allocate solver dof to variable maps equations indices.", &
                                    & err,error,*999)
                                  ALLOCATE(solverMapping%solverColToEquationsColsMap(solverMatrixIdx)% &
                                    & solverDofToVariableMaps(solverLocalDof(rank))%variable(numberColEquationsCols),STAT=err)
                                  IF(err/=0) CALL FlagError("Could not allocate solver dof to variable maps variable type.", &
                                    & err,error,*999)
                                  ALLOCATE(solverMapping%solverColToEquationsColsMap(solverMatrixIdx)% &
                                    & solverDofToVariableMaps(solverLocalDof(rank))%variableDof(numberColEquationsCols),STAT=err)
                                  IF(err/=0) CALL FlagError("Could not allocate solver dof to variable maps variable dof.", &
                                    & err,error,*999)
                                  ALLOCATE(solverMapping%solverColToEquationsColsMap(solverMatrixIdx)% &
                                    & solverDofToVariableMaps(solverLocalDof(rank))%variableCoefficient(numberColEquationsCols), &
                                    & STAT=err)
                                  IF(err/=0) &
                                    & CALL FlagError("Could not allocate solver dof to variable maps variable coefficient.", &
                                    & err,error,*999)
                                  ALLOCATE(solverMapping%solverColToEquationsColsMap(solverMatrixIdx)% &
                                    & solverDofToVariableMaps(solverLocalDof(rank))%additiveConstant(numberColEquationsCols), &
                                    & STAT=err)
                                  IF(err/=0) CALL FlagError("Could not allocate solver dof to variable maps additive constant.", &
                                    & err,error,*999)
                                  !Setup
                                  solverMapping%solverColToEquationsColsMap(solverMatrixIdx)%solverDofToVariableMaps( &
                                    & solverLocalDof(rank))%numberOfEquations=numberColEquationsCols
                                  DO colEquationsColIdx=1,numberColEquationsCols
                                    eqnLocalDof=colEquationCols%localDofs(colEquationsColIdx)
                                    couplingCoefficient=colEquationCols%coefficients(colEquationsColIdx)
                                    solverMapping%solverColToEquationsColsMap(solverMatrixIdx)%solverDofToVariableMaps( &
                                      & solverLocalDof(rank))%equationsTypes(colEquationsColIdx)= &
                                      & SOLVER_MAPPING_EQUATIONS_INTERFACE_CONDITION
                                    solverMapping%solverColToEquationsColsMap(solverMatrixIdx)%solverDofToVariableMaps( &
                                      & solverLocalDof(rank))%equationsIndices(colEquationsColIdx)=interfaceConditionIdx
                                    solverMapping%solverColToEquationsColsMap(solverMatrixIdx)%solverDofToVariableMaps( &
                                      & solverLocalDof(rank))%variable(colEquationsColIdx)%ptr=>lagrangeVariable
                                    solverMapping%solverColToEquationsColsMap(solverMatrixIdx)%solverDofToVariableMaps( &
                                      & solverLocalDof(rank))%variableDof(colEquationsColIdx)=eqnLocalDof
                                    solverMapping%solverColToEquationsColsMap(solverMatrixIdx)%solverDofToVariableMaps( &
                                      & solverLocalDof(rank))%variableCoefficient(colEquationsColIdx)=couplingCoefficient
                                    solverMapping%solverColToEquationsColsMap(solverMatrixIdx)%solverDofToVariableMaps( &
                                      & solverLocalDof(rank))%additiveConstant(colEquationsColIdx)=0.0_DP
                                  ENDDO !colEquationsColIdx
                                ENDIF
                                DO colEquationsColIdx=1,numberColEquationsCols
                                  eqnLocalDof=colEquationCols%localDofs(colEquationsColIdx)
                                  couplingCoefficient=colEquationCols%coefficients(colEquationsColIdx)
                                  IF(equationType==SOLVER_MAPPING_EQUATIONS_INTERFACE_CONDITION) THEN
                                    IF(.NOT.variableRankProcessed(solverVariableIdx,rank)) THEN
                                      !Set up the equations variables -> solver columns mapping
                                      solverMapping%interfaceConditionToSolverMap(interfaceConditionIdx)% &
                                        & interfaceToSolverMatrixMapsSm(solverMatrixIdx)%lagrangeVariableToSolverColMap% &
                                        & columnNumbers(eqnLocalDof)=solverGlobalDof
                                      solverMapping%interfaceConditionToSolverMap(interfaceConditionIdx)% &
                                        & interfaceToSolverMatrixMapsSm(solverMatrixIdx)%lagrangeVariableToSolverColMap% &
                                        & couplingCoefficients(eqnLocalDof)=couplingCoefficient
                                      solverMapping%interfaceConditionToSolverMap(interfaceConditionIdx)% &
                                        & interfaceToSolverMatrixMapsSm(solverMatrixIdx)%lagrangeVariableToSolverColMap% &
                                        & additiveConstants(eqnLocalDof)=0.0_DP
                                      !Set up the equations columns -> solver columns mapping
                                      interfaceColumn=interfaceMapping%LAGRANGE_DOF_TO_COLUMN_MAP(localDof)
                                      !Allocate the equation to solver map column items.
                                      ALLOCATE(solverMapping%interfaceConditionToSolverMap(interfaceConditionIdx)% &
                                        & interfaceToSolverMatrixMapsSm(solverMatrixIdx)%interfaceColToSolverColsMap( &
                                        & interfaceColumn)%solverCols(1),STAT=err)
                                      IF(err/=0) CALL  FlagError("Could not allocate interface column to solver columns map "// &
                                        & "solver columns.",err,error,*999)
                                      ALLOCATE(solverMapping%interfaceConditionToSolverMap(interfaceConditionIdx)% &
                                        & interfaceToSolverMatrixMapsSm(solverMatrixIdx)%interfaceColToSolverColsMap( &
                                        & interfaceColumn)%couplingCoefficients(1),STAT=err)
                                      IF(err/=0) CALL FlagError("Could not allocate interface column to solver columns map "// &
                                        & "coupling coefficients.",err,error,*999)
                                      solverMapping%interfaceConditionToSolverMap(interfaceConditionIdx)% &
                                        & interfaceToSolverMatrixMapsSm(solverMatrixIdx)%interfaceColToSolverColsMap( &
                                        & interfaceColumn)%numberOfSolverCols=1
                                      solverMapping%interfaceConditionToSolverMap(interfaceConditionIdx)% &
                                        & interfaceToSolverMatrixMapsSm(solverMatrixIdx)%interfaceColToSolverColsMap( &
                                        & interfaceColumn)%solverCols(1)=solverGlobalDof
                                      solverMapping%interfaceConditionToSolverMap(interfaceConditionIdx)% &
                                        & interfaceToSolverMatrixMapsSm(solverMatrixIdx)%interfaceColToSolverColsMap( &
                                        & interfaceColumn)%couplingCoefficients(1)=1.0_DP
                                    ENDIF
                                  ELSE
                                    !Set up the equations variables -> solver columns mapping
                                    solverMapping%interfaceConditionToSolverMap(interfaceConditionIdx)% &
                                      & interfaceToSolverMatrixMapsSm(solverMatrixIdx)%dependentVariableToSolverColMaps( &
                                      & interfaceMatrixIdx)%columnNumbers(eqnLocalDof)=solverGlobalDof
                                    solverMapping%interfaceConditionToSolverMap(interfaceConditionIdx)% &
                                      & interfaceToSolverMatrixMapsSm(solverMatrixIdx)%dependentVariableToSolverColMaps( &
                                      & interfaceMatrixIdx)%couplingCoefficients(eqnLocalDof)=1.0_DP
                                    solverMapping%interfaceConditionToSolverMap(interfaceConditionIdx)% &
                                      & interfaceToSolverMatrixMapsSm(solverMatrixIdx)%dependentVariableToSolverColMaps( &
                                      & interfaceMatrixIdx)%additiveConstants(eqnLocalDof)=0.0_DP
                                    !Set up the equations columns -> solver columns mapping
                                    interfaceRow=interfaceMapping%INTERFACE_MATRIX_ROWS_TO_VAR_MAPS(interfaceMatrixIdx)% &
                                      & VARIABLE_DOF_TO_ROW_MAP(eqnLocalDof)
                                    !Allocate the equation to solver map column items.
                                    ALLOCATE(solverMapping%interfaceConditionToSolverMap(interfaceConditionIdx)% &
                                      & interfaceToSolverMatrixMapsSm(solverMatrixIdx)% &
                                      & interfaceEquationsToSolverMatrixMaps(interfaceMatrixIdx)%ptr% &
                                      & interfaceRowToSolverColsMap(interfaceRow)%solverCols(1),STAT=err)
                                    IF(err/=0) &
                                      & CALL FlagError("Could not allocate interface equations row to solver columns map "// &
                                      & "solver columns.",err,error,*999)
                                    ALLOCATE(solverMapping%interfaceConditionToSolverMap(interfaceConditionIdx)% &
                                      & interfaceToSolverMatrixMapsSm(solverMatrixIdx)% &
                                      & interfaceEquationsToSolverMatrixMaps(interfaceMatrixIdx)%ptr% &
                                      & interfaceRowToSolverColsMap(interfaceRow)%couplingCoefficients(1),STAT=err)
                                    IF(err/=0) &
                                      & CALL FlagError("Could not allocate interface equations row to solver columns map "// &
                                      & "coupling coefficients.",err,error,*999)
                                    solverMapping%interfaceConditionToSolverMap(interfaceConditionIdx)% &
                                      & interfaceToSolverMatrixMapsSm(solverMatrixIdx)% &
                                      & interfaceEquationsToSolverMatrixMaps(interfaceMatrixIdx)%ptr% &
                                      & interfaceRowToSolverColsMap(interfaceRow)%numberOfSolverCols=1
                                    solverMapping%interfaceConditionToSolverMap(interfaceConditionIdx)% &
                                      & interfaceToSolverMatrixMapsSm(solverMatrixIdx)% &
                                      & interfaceEquationsToSolverMatrixMaps(interfaceMatrixIdx)%ptr% &
                                      & interfaceRowToSolverColsMap(interfaceRow)%solverCols(1)=solverGlobalDof
                                    solverMapping%interfaceConditionToSolverMap(interfaceConditionIdx)% &
                                      & interfaceToSolverMatrixMapsSm(solverMatrixIdx)% &
                                      & interfaceEquationsToSolverMatrixMaps(interfaceMatrixIdx)%ptr% &
                                      & interfaceRowToSolverColsMap(interfaceRow)%couplingCoefficients(1)= &
                                      & interfaceMapping%INTERFACE_MATRIX_ROWS_TO_VAR_MAPS(interfaceMatrixIdx)%MATRIX_COEFFICIENT
                                  ENDIF
                                ENDDO !colEquationsColIdx
                              ELSE !rank /= myrank
                                IF(.NOT.variableRankProcessed(solverVariableIdx,rank)) THEN
                                  
                                  !Set up the column domain mappings.
                                  CALL DOMAIN_MAPPINGS_MAPPING_GLOBAL_INITIALISE(colDomainMapping%GLOBAL_TO_LOCAL_MAP( &
                                    & solverGlobalDof),err,error,*999)
                                  !There are no ghosted cols for the solver matrices so there is only 1 domain for the global to
                                  !local map.
                                  !Allocate the global to local map arrays
                                  ALLOCATE(colDomainMapping%GLOBAL_TO_LOCAL_MAP(solverGlobalDof)%LOCAL_NUMBER(1),STAT=err)
                                  IF(err/=0) &
                                    & CALL FlagError("Could not allocate column domain global to local map local number.", &
                                    & err,error,*999)
                                  ALLOCATE(colDomainMapping%GLOBAL_TO_LOCAL_MAP(solverGlobalDof)%DOMAIN_NUMBER(1),STAT=err)
                                  IF(err/=0) &
                                    & CALL FlagError("Could not allocate column domain global to local map domain number.", &
                                    & err,error,*999)
                                  ALLOCATE(colDomainMapping%GLOBAL_TO_LOCAL_MAP(solverGlobalDof)%LOCAL_TYPE(1),STAT=err)
                                  IF(err/=0) &
                                    & CALL FlagError("Could not allocate column domain global to local map domain number.", &
                                    & err,error,*999)
                                  !Set up the global to local mappings
                                  colDomainMapping%GLOBAL_TO_LOCAL_MAP(solverGlobalDof)%NUMBER_OF_DOMAINS=1
                                  colDomainMapping%GLOBAL_TO_LOCAL_MAP(solverGlobalDof)%LOCAL_NUMBER(1)=solverLocalDof(rank)
                                  colDomainMapping%GLOBAL_TO_LOCAL_MAP(solverGlobalDof)%DOMAIN_NUMBER(1)=rank
                                  colDomainMapping%GLOBAL_TO_LOCAL_MAP(solverGlobalDof)%LOCAL_TYPE(1)=DOMAIN_LOCAL_INTERNAL
                                  
                                ENDIF
                              ENDIF !rank == myrank
                            ELSE IF(constrainedDof) THEN
                              !Do nothing, this is set up above
                            ELSE
                              IF(rank==myrank) THEN
                                IF(equationType==SOLVER_MAPPING_EQUATIONS_INTERFACE_CONDITION) THEN
                                  !Set up the equations variables -> solver columns mapping
                                  solverMapping%interfaceConditionToSolverMap(interfaceConditionIdx)% &
                                    & interfaceToSolverMatrixMapsSm(solverMatrixIdx)% &
                                    & lagrangeVariableToSolverColMap%columnNumbers(localDof)=0
                                  solverMapping%interfaceConditionToSolverMap(interfaceConditionIdx)% &
                                    & interfaceToSolverMatrixMapsSm(solverMatrixIdx)% &
                                    & lagrangeVariableToSolverColMap%couplingCoefficients(localDof)=0.0_DP
                                  solverMapping%interfaceConditionToSolverMap(interfaceConditionIdx)% &
                                    & interfaceToSolverMatrixMapsSm(solverMatrixIdx)% &
                                    & lagrangeVariableToSolverColMap%additiveConstants(localDof)=0.0_DP
                                  interfaceColumn=interfaceMapping%LAGRANGE_DOF_TO_COLUMN_MAP(localDof)
                                  solverMapping%interfaceConditionToSolverMap(interfaceConditionIdx)% &
                                    & interfaceToSolverMatrixMapsSm(solverMatrixIdx)%interfaceColToSolverColsMap( &
                                    & interfaceColumn)%numberOfSolverCols=0
                                ELSE
                                  !Set up the equations variables -> solver columns mapping
                                  solverMapping%interfaceConditionToSolverMap(interfaceConditionIdx)% &
                                    & interfaceToSolverMatrixMapsSm(solverMatrixIdx)% &
                                    & dependentVariableToSolverColMaps(interfaceMatrixIdx)%columnNumbers(localDof)=0
                                  solverMapping%interfaceConditionToSolverMap(interfaceConditionIdx)% &
                                    & interfaceToSolverMatrixMapsSm(solverMatrixIdx)% &
                                    & dependentVariableToSolverColMaps(interfaceMatrixIdx)% &
                                    & couplingCoefficients(localDof)=0.0_DP
                                  solverMapping%interfaceConditionToSolverMap(interfaceConditionIdx)% &
                                    & interfaceToSolverMatrixMapsSm(solverMatrixIdx)% &
                                    & dependentVariableToSolverColMaps(interfaceMatrixIdx)% &
                                    & additiveConstants(localDof)=0.0_DP
                                  !Set up the equations columns -> solver columns mapping
                                  interfaceRow=interfaceMapping%INTERFACE_MATRIX_ROWS_TO_VAR_MAPS(interfaceMatrixIdx)% &
                                    & VARIABLE_DOF_TO_ROW_MAP(localDof)
                                  solverMapping%interfaceConditionToSolverMap(interfaceConditionIdx)% &
                                    & interfaceToSolverMatrixMapsSm(solverMatrixIdx)% &
                                    & interfaceEquationsToSolverMatrixMaps(interfaceMatrixIdx)%ptr% &
                                    & interfaceRowToSolverColsMap(interfaceRow)%numberOfSolverCols=0
                                ENDIF
                              ENDIF !rank==myrank
                            ENDIF !include_column
                          ENDDO !globalDof
                        CASE(INTERFACE_CONDITION_AUGMENTED_LAGRANGE_METHOD)
                          CALL FlagError("Not implemented.",err,error,*999)
                        CASE(INTERFACE_CONDITION_POINT_TO_POINT_METHOD)
                          CALL FlagError("Not implemented.",err,error,*999)
                        CASE DEFAULT
                          localError="The interface condition method of "// &
                            & TRIM(NumberToVString(interfaceCondition%METHOD,"*",err,error))// &
                            & " is invalid."
                          CALL FlagError(localError,err,error,*999)
                        END SELECT
                      CASE DEFAULT
                        localError="The equation type of "//TRIM(NumberToVString(equationType,"*",err,error))// &
                          & " is invalid."
                        CALL FlagError(localError,err,error,*999)
                      END SELECT
                      variableRankProcessed(solverVariableIdx,rank)=.TRUE.
                    ENDIF !Number of rank columns > 0
                    IF(ALLOCATED(rankGlobalColsList)) DEALLOCATE(rankGlobalColsList)              
                  ENDDO !equation_idx
                  
                ENDDO !solverVariableIdx
              ENDDO !rank
            ENDDO !dofType
            
            !Deallocate dof map
            DO solverVariableIdx=1,solverMapping%columnVariablesList(solverMatrixIdx)%numberOfVariables
              DEALLOCATE(dofMap(solverVariableIdx)%ptr)
            ENDDO !solverVariableIdx
            DEALLOCATE(dofMap)
            !Deallocate solver local dof
            IF(ALLOCATED(solverLocalDof)) DEALLOCATE(solverLocalDof)
            
            CALL DOMAIN_MAPPINGS_LOCAL_FROM_GLOBAL_CALCULATE(colDomainMapping,err,error,*999)
            
           IF(ALLOCATED(subMatrixInformation)) DEALLOCATE(subMatrixInformation)
            IF(ALLOCATED(subMatrixList)) DEALLOCATE(subMatrixList)
            IF(ALLOCATED(variableRankProcessed)) DEALLOCATE(variableRankProcessed)
            IF(ALLOCATED(numberOfVariableGlobalSolverDofs)) DEALLOCATE(numberOfVariableGlobalSolverDofs)
            IF(ALLOCATED(numberOfVariableLocalSolverDofs)) DEALLOCATE(numberOfVariableLocalSolverDofs)
            IF(ALLOCATED(totalNumberOfVariableLocalSolverDofs)) DEALLOCATE(totalNumberOfVariableLocalSolverDofs)
            DO solverVariableIdx=1,numberOfEquationsVariables+numberOfInterfaceVariables
              CALL LIST_DESTROY(variablesList(solverVariableIdx)%ptr,err,error,*999)
            ENDDO !solverVariableIdx
          ENDDO !solverMatrixIdx
          IF(ALLOCATED(dummyDofCoupling%globalDofs)) DEALLOCATE(dummyDofCoupling%globalDofs)
          IF(ALLOCATED(dummyDofCoupling%localDofs)) DEALLOCATE(dummyDofCoupling%localDofs)
          IF(ALLOCATED(dummyDofCoupling%coefficients)) DEALLOCATE(dummyDofCoupling%coefficients)
          CALL SolverDofCouplings_Finalise(columnCouplings,err,error,*999)

          !
          ! 5. Set up the column mappings such that the solver matrix (sm) and equations matrix (em) orderings are the same.
          !
          
          DO equationsSetIdx=1,solverMapping%numberOfEquationsSets
            equationsSet=>solverMapping%equationsSets(equationsSetIdx)%ptr
            equations=>equationsSet%equations
            equationsMapping=>equations%EQUATIONS_MAPPING
            dynamicMapping=>equationsMapping%DYNAMIC_MAPPING
            linearMapping=>equationsMapping%LINEAR_MAPPING
            IF(ASSOCIATED(dynamicMapping)) THEN
              DO equationsMatrixIdx=1,dynamicMapping%NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES
                ALLOCATE(solverMapping%equationsSetToSolverMap(equationsSetIdx)%equationsToSolverMatrixMapsEm( &
                  & equationsMatrixIdx)%equationsToSolverMatrixMaps(solverMapping%equationsSetToSolverMap( &
                  & equationsSetIdx)%equationsToSolverMatrixMapsEm(equationsMatrixIdx)%numberOfSolverMatrices), &
                  & STAT=err)
                IF(err/=0) CALL FlagError("Could not allocate equations to solver matrix maps.",err,error,*999)
                solverMapping%equationsSetToSolverMap(equationsSetIdx)%equationsToSolverMatrixMapsEm( &
                  & equationsMatrixIdx)%numberOfSolverMatrices=0
              ENDDO !equationsMatrixIdx
              DO solverMatrixIdx=1,solverMapping%numberOfSolverMatrices
                DO equationsMatrixIdx=1,solverMapping%equationsSetToSolverMap(equationsSetIdx)% &
                  & equationsToSolverMatrixMapsSm(solverMatrixIdx)%numberOfDynamicEquationsMatrices
                  IF(solverMapping%equationsSetToSolverMap(equationsSetIdx)%equationsToSolverMatrixMapsSm( &
                    & solverMatrixIdx)%dynamicEquationsToSolverMatrixMaps(equationsMatrixIdx)%ptr% &
                    & solverMatrixNumber==solverMatrixIdx) THEN
                    solverMapping%equationsSetToSolverMap(equationsSetIdx)%equationsToSolverMatrixMapsEm( &
                      & equationsMatrixIdx)%numberOfSolverMatrices=solverMapping%equationsSetToSolverMap( &
                      & equationsSetIdx)%equationsToSolverMatrixMapsEm(equationsMatrixIdx)%numberOfSolverMatrices+1
                    solverMapping%equationsSetToSolverMap(equationsSetIdx)%equationsToSolverMatrixMapsEm( &
                      & equationsMatrixIdx)%equationsToSolverMatrixMaps(solverMapping%equationsSetToSolverMap( &
                      & equationsSetIdx)%equationsToSolverMatrixMapsEm(equationsMatrixIdx)%numberOfSolverMatrices)% &
                      & PTR=>solverMapping%equationsSetToSolverMap(equationsSetIdx)%equationsToSolverMatrixMapsSm( &
                      & solverMatrixIdx)%dynamicEquationsToSolverMatrixMaps(equationsMatrixIdx)%ptr
                  ENDIF
                ENDDO !equationsMatrixIdx
              ENDDO !solverMatrixIdx             
            ELSE IF(ASSOCIATED(linearMapping)) THEN
              DO equationsMatrixIdx=1,linearMapping%NUMBER_OF_LINEAR_EQUATIONS_MATRICES
                ALLOCATE(solverMapping%equationsSetToSolverMap(equationsSetIdx)%equationsToSolverMatrixMapsEm( &
                  & equationsMatrixIdx)%equationsToSolverMatrixMaps(solverMapping%equationsSetToSolverMap( &
                  & equationsSetIdx)%equationsToSolverMatrixMapsEm(equationsMatrixIdx)%numberOfSolverMatrices), &
                  & STAT=err)
                IF(err/=0) CALL FlagError("Could not allocate equations to solver matrix maps.",err,error,*999)
                solverMapping%equationsSetToSolverMap(equationsSetIdx)%equationsToSolverMatrixMapsEm( &
                  & equationsMatrixIdx)%numberOfSolverMatrices=0
              ENDDO !equationsMatrixIdx
              DO solverMatrixIdx=1,solverMapping%numberOfSolverMatrices
                DO equationsMatrixIdx=1,solverMapping%equationsSetToSolverMap(equationsSetIdx)% &
                  & equationsToSolverMatrixMapsSm(solverMatrixIdx)%numberOfLinearEquationsMatrices
                  IF(solverMapping%equationsSetToSolverMap(equationsSetIdx)%equationsToSolverMatrixMapsSm( &
                    & solverMatrixIdx)%linearEquationsToSolverMatrixMaps(equationsMatrixIdx)%ptr% &
                    & solverMatrixNumber==solverMatrixIdx) THEN
                    solverMapping%equationsSetToSolverMap(equationsSetIdx)%equationsToSolverMatrixMapsEm( &
                      & equationsMatrixIdx)%numberOfSolverMatrices=solverMapping%equationsSetToSolverMap( &
                      & equationsSetIdx)%equationsToSolverMatrixMapsEm(equationsMatrixIdx)%numberOfSolverMatrices+1
                    solverMapping%equationsSetToSolverMap(equationsSetIdx)%equationsToSolverMatrixMapsEm( &
                      & equationsMatrixIdx)%equationsToSolverMatrixMaps(solverMapping%equationsSetToSolverMap( &
                      & equationsSetIdx)%equationsToSolverMatrixMapsEm(equationsMatrixIdx)%numberOfSolverMatrices)% &
                      & PTR=>solverMapping%equationsSetToSolverMap(equationsSetIdx)%equationsToSolverMatrixMapsSm( &
                      & solverMatrixIdx)%linearEquationsToSolverMatrixMaps(equationsMatrixIdx)%ptr
                  ENDIF
                ENDDO !equationsMatrixIdx
              ENDDO !solverMatrixIdx             
            ENDIF
          ENDDO !equationsSetIdx
          DO interfaceConditionIdx=1,solverMapping%numberOfInterfaceConditions
            interfaceCondition=>solverMapping%interfaceConditions(interfaceConditionIdx)%ptr
            SELECT CASE(interfaceCondition%method)
            CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD)
              interfaceEquations=>interfaceCondition%INTERFACE_EQUATIONS
              interfaceMapping=>interfaceEquations%INTERFACE_MAPPING
              DO interfaceMatrixIdx=1,interfaceMapping%NUMBER_OF_INTERFACE_MATRICES
                ALLOCATE(solverMapping%interfaceConditionToSolverMap(interfaceConditionIdx)% &
                  & interfaceToSolverMatrixMapsIm(interfaceMatrixIdx)%interfaceToSolverMatrixMaps(solverMapping% &
                  & interfaceConditionToSolverMap(interfaceConditionIdx)%interfaceToSolverMatrixMapsIm( &
                  & interfaceMatrixIdx)%numberOfSolverMatrices),STAT=err)
                IF(err/=0) CALL FlagError("Could not allocate interface to solver matrix maps.",err,error,*999)
                solverMapping%interfaceConditionToSolverMap(interfaceConditionIdx)%interfaceToSolverMatrixMapsIm( &
                  & interfaceMatrixIdx)%numberOfSolverMatrices=0
              ENDDO !interfaceMatrixIdx
              DO solverMatrixIdx=1,solverMapping%numberOfSolverMatrices
                DO interfaceMatrixIdx=1,solverMapping%interfaceConditionToSolverMap(interfaceConditionIdx)% &
                  & interfaceToSolverMatrixMapsSm(solverMatrixIdx)%numberOfInterfaceMatrices
                  IF(solverMapping%interfaceConditionToSolverMap(interfaceConditionIdx)%interfaceToSolverMatrixMapsSm( &
                    & solverMatrixIdx)%interfaceEquationsToSolverMatrixMaps(interfaceMatrixIdx)%ptr% &
                    & solverMatrixNumber==solverMatrixIdx) THEN
                    solverMapping%interfaceConditionToSolverMap(interfaceConditionIdx)%interfaceToSolverMatrixMapsIm( &
                      & interfaceMatrixIdx)%numberOfSolverMatrices=solverMapping%interfaceConditionToSolverMap( &
                      & interfaceConditionIdx)%interfaceToSolverMatrixMapsIm(interfaceMatrixIdx)% &
                      & numberOfSolverMatrices+1
                    solverMapping%interfaceConditionToSolverMap(interfaceConditionIdx)%interfaceToSolverMatrixMapsIm( &
                      & interfaceMatrixIdx)%interfaceToSolverMatrixMaps(solverMapping%interfaceConditionToSolverMap( &
                      & interfaceConditionIdx)%interfaceToSolverMatrixMapsIm(interfaceMatrixIdx)%&
                      & numberOfSolverMatrices)%ptr=>solverMapping%interfaceConditionToSolverMap(interfaceConditionIdx)% &
                      & interfaceToSolverMatrixMapsSm(solverMatrixIdx)%interfaceEquationsToSolverMatrixMaps( &
                      & interfaceMatrixIdx)%ptr
                  ENDIF
                ENDDO !interfaceMatrixIdx
              ENDDO !solverMatrixIdx             
            CASE(INTERFACE_CONDITION_AUGMENTED_LAGRANGE_METHOD)
              CALL FlagError("Not implemented.",err,error,*999)
            CASE(INTERFACE_CONDITION_POINT_TO_POINT_METHOD)
              CALL FlagError("Not implemented.",err,error,*999)
            CASE DEFAULT
              localError="The interface condition method of "// &
                & TRIM(NumberToVString(interfaceCondition%METHOD,"*",err,error))//" is invalid."
              CALL FlagError(localError,err,error,*999)
            END SELECT
                        
          ENDDO !interfaceConditionIdx
        ELSE
          CALL FlagError("The solver mapping solver equations are not associated.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Solver mapping create values cache is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Solver mapping is not associated.",err,error,*999)
    ENDIF

    IF(DIAGNOSTICS1) THEN
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"Solver mappings:",err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Number of solver matrices = ",solverMapping%numberOfSolverMatrices, &
        & err,error,*999)               
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Equation sets:",err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Number of equations sets = ",solverMapping%numberOfEquationsSets, &
        & err,error,*999)
      DO equationsSetIdx=1,solverMapping%numberOfEquationsSets
        equationsSet=>solverMapping%equationsSets(equationsSetIdx)%ptr
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Equations set index : ",equationsSetIdx,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Region user number        = ",equationsSet%region%USER_NUMBER, &
          & err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Equations set user number = ",equationsSet%USER_NUMBER, &
          & err,error,*999)                
      ENDDO !equationsSetIdx
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Interface conditions:",err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Number of interface conditions = ",solverMapping% &
        & numberOfInterfaceConditions,err,error,*999)
      DO interfaceConditionIdx=1,solverMapping%numberOfInterfaceConditions
        interfaceCondition=>solverMapping%interfaceConditions(interfaceConditionIdx)%ptr
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Interface condition index : ",interfaceConditionIdx, &
          & err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Parent region user number       = ",interfaceCondition%interface% &
          & PARENT_REGION%USER_NUMBER,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Interface condition user number = ",interfaceCondition%USER_NUMBER, &
          & err,error,*999)                
      ENDDO !equationsSetIdx
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Equations variables list:",err,error,*999)
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"    Row variables list:",err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Number of variables = ",solverMapping%rowVariablesList% &
        & numberOfVariables,err,error,*999)        
      DO variableIdx=1,solverMapping%rowVariablesList%numberOfVariables
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Variable : ",variableIdx,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Variable type = ",solverMapping%rowVariablesList% &
          & variables(variableIdx)%variableType,err,error,*999)        
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Number of equations = ",solverMapping%rowVariablesList% &
          & variables(variableIdx)%numberOfEquations,err,error,*999)        
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,solverMapping%rowVariablesList%variables(variableIdx)% &
          & numberOfEquations,5,5,solverMapping%rowVariablesList%variables(variableIdx)%equationTypes, &
          & '("        Equation types   :",5(X,I13))','(26X,5(X,I13))',err,error,*999)
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,solverMapping%rowVariablesList%variables(variableIdx)% &
          & numberOfEquations,5,5,solverMapping%rowVariablesList%variables(variableIdx)%equationTypes, &
          & '("        Equation indices :",5(X,I13))','(26X,5(X,I13))',err,error,*999)
      ENDDO !variableIdx
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"    Column variables lists:",err,error,*999)
      DO solverMatrixIdx=1,solverMapping%numberOfSolverMatrices
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Solver matrix : ",solverMatrixIdx,err,error,*999)        
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Number of variables = ",solverMapping%columnVariablesList( &
          & solverMatrixIdx)%numberOfVariables,err,error,*999)        
        DO variableIdx=1,solverMapping%columnVariablesList(solverMatrixIdx)%numberOfVariables
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Variable : ",variableIdx,err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"          Variable type = ",solverMapping%columnVariablesList( &
            & solverMatrixIdx)%variables(variableIdx)%variableType,err,error,*999)        
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"          Number of equations = ",solverMapping%columnVariablesList( &
            & solverMatrixIdx)%variables(variableIdx)%numberOfEquations,err,error,*999)        
          CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,solverMapping%columnVariablesList(solverMatrixIdx)%variables( &
            & variableIdx)%numberOfEquations,5,5,solverMapping%columnVariablesList(solverMatrixIdx)%variables( &
            & variableIdx)%equationTypes,'("          Equation types   :",5(X,I13))','(28X,5(X,I13))',err,error,*999)
          CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,solverMapping%columnVariablesList(solverMatrixIdx)%variables( &
            & variableIdx)%numberOfEquations,5,5,solverMapping%columnVariablesList(solverMatrixIdx)%variables( &
            & variableIdx)%equationTypes,'("          Equation indices :",5(X,I13))','(28X,5(X,I13))',err,error,*999)
        ENDDO !variableIdx
      ENDDO !solverMatrixIdx
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Solver row to equations rows mappings:",err,error,*999)      
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Number of rows = ",solverMapping%numberOfRows,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Number of global rows = ",solverMapping%numberOfGlobalRows, &
        & err,error,*999)
      IF(DIAGNOSTICS2) THEN
        DO rowIdx=1,solverMapping%numberOfRows
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Solver row : ",rowIdx,err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Number of equations sets mapped to = ",solverMapping% &
            & solverRowToEquationsRowsMap(rowIdx)%numberOfEquationsSets,err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Interface condition index          = ",solverMapping% &
            & solverRowToEquationsRowsMap(rowIdx)%interfaceConditionIndex,err,error,*999)
          IF(solverMapping%solverRowToEquationsRowsMap(rowIdx)%interfaceConditionIndex==0) THEN
            !Row is an equations set row
            CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,solverMapping%solverRowToEquationsRowsMap(rowIdx)% &
              & numberOfEquationsSets,5,5,solverMapping%solverRowToEquationsRowsMap(rowIdx)%equationsIndex, &
              & '("      Equations indices      :",5(X,I13))','(30X,5(X,I13))',err,error,*999) 
            CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,solverMapping%solverRowToEquationsRowsMap(rowIdx)% &
              & numberOfEquationsSets,5,5,solverMapping%solverRowToEquationsRowsMap(rowIdx)%rowcolNumber, &
              & '("      Equations row numbers  :",5(X,I13))','(30X,5(X,I13))',err,error,*999) 
            CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,solverMapping%solverRowToEquationsRowsMap(rowIdx)% &
              & numberOfEquationsSets,5,5,solverMapping%solverRowToEquationsRowsMap(rowIdx)%couplingCoefficients, &
              & '("      Coupling coefficients  :",5(X,E13.6))','(30X,5(X,E13.6))',err,error,*999)
          ELSE
            !Row is an interface condition row
!!TODO: format better
            CALL WRITE_STRING_FMT_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Interface col numbers : ",solverMapping% &
              & solverRowToEquationsRowsMap(rowIdx)%rowcolNumber(1),"(I13)",err,error,*999) 
            CALL WRITE_STRING_FMT_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Coupling coefficients : ",solverMapping% &
              & solverRowToEquationsRowsMap(rowIdx)%couplingCoefficients(1),"(E13.6)",err,error,*999)
          ENDIF
        ENDDO !rowIdx
      ENDIF
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Solver column to equations column mappings:",err,error,*999)      
      DO solverMatrixIdx=1,solverMapping%numberOfSolverMatrices
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Solver matrix : ",solverMatrixIdx,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Number of columns = ",solverMapping% &
          & solverColToEquationsColsMap(solverMatrixIdx)%numberOfColumns,err,error,*999)
        IF(DIAGNOSTICS2) THEN
          CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"      Solver column to equations set columns mappings:",err,error,*999)
          DO equationsSetIdx=1,solverMapping%numberOfEquationsSets
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Equations set index : ",equationsSetIdx,err,error,*999)
            DO columnIdx=1,solverMapping%solverColToEquationsColsMap(solverMatrixIdx)%numberOfColumns           
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE, "          Solver column : ",columnIdx,err,error,*999)
              IF(solverMapping%solverColToEquationsColsMap(solverMatrixIdx)%solverColToEquationsSetMaps( &
                & equationsSetIdx)%haveDynamic) THEN
               CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"            Number of dynamic equations matrices mapped to = ", &
                  & solverMapping%solverColToEquationsColsMap(solverMatrixIdx)%solverColToEquationsSetMaps( &
                  & equationsSetIdx)%solverColToDynamicEquationsMaps(columnIdx)%numberOfDynamicEquationsMatrices, &
                  & err,error,*999)
                IF(solverMapping%solverColToEquationsColsMap(solverMatrixIdx)%solverColToEquationsSetMaps( &
                  & equationsSetIdx)%solverColToDynamicEquationsMaps(columnIdx)%numberOfDynamicEquationsMatrices>0) THEN
                  CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,solverMapping%solverColToEquationsColsMap( &
                    & solverMatrixIdx)%solverColToEquationsSetMaps(equationsSetIdx)%solverColToDynamicEquationsMaps( &
                    & columnIdx)%numberOfDynamicEquationsMatrices,5,5,solverMapping%solverColToEquationsColsMap( &
                    & solverMatrixIdx)%solverColToEquationsSetMaps(equationsSetIdx)%solverColToDynamicEquationsMaps( &
                    & columnIdx)%equationsMatrixNumbers,'("            Equation matrices numbers :",5(X,I13))', &
                    & '(39X,5(X,I13))',err,error,*999)
                  CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,solverMapping%solverColToEquationsColsMap( &
                    & solverMatrixIdx)%solverColToEquationsSetMaps(equationsSetIdx)%solverColToDynamicEquationsMaps( &
                    & columnIdx)%numberOfDynamicEquationsMatrices,5,5,solverMapping%solverColToEquationsColsMap( &
                    & solverMatrixIdx)%solverColToEquationsSetMaps(equationsSetIdx)%solverColToDynamicEquationsMaps( &
                    & columnIdx)%equationsColNumbers,'("            Equation column numbers   :",5(X,I13))', &
                    & '(39X,5(X,I13))',err,error,*999)
                  CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,solverMapping%solverColToEquationsColsMap( &
                    & solverMatrixIdx)%solverColToEquationsSetMaps(equationsSetIdx)%solverColToDynamicEquationsMaps( &
                    & columnIdx)%numberOfDynamicEquationsMatrices,5,5,solverMapping%solverColToEquationsColsMap( &
                    & solverMatrixIdx)%solverColToEquationsSetMaps(equationsSetIdx)%solverColToDynamicEquationsMaps( &
                    & columnIdx)%couplingCoefficients,'("            Coupling coefficients     :",5(X,E13.6))', &
                    & '(39X,5(X,E13.6))',err,error,*999)
                ENDIF
              ELSE
                 CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"            Number of dynamic equations matrices mapped to = ", &
                  & 0_INTG,err,error,*999)
              ENDIF
!!TODO what about dynamic nonlinear mappings???
              IF(solverMapping%solverColToEquationsColsMap(solverMatrixIdx)%solverColToEquationsSetMaps( &
                & equationsSetIdx)%haveStatic) THEN
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE, "            Number of linear equations matrices mapped to  = ", &
                  & solverMapping%solverColToEquationsColsMap(solverMatrixIdx)%solverColToEquationsSetMaps( &
                  & equationsSetIdx)%SolverColToStaticEquationsMaps(columnIdx)%numberOfLinearEquationsMatrices, &
                  & err,error,*999)
                IF(solverMapping%solverColToEquationsColsMap(solverMatrixIdx)%solverColToEquationsSetMaps( &
                  & equationsSetIdx)%SolverColToStaticEquationsMaps(columnIdx)%numberOfLinearEquationsMatrices>0) THEN
                  !CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,solverMapping%solverColToEquationsColsMap( &
                  !  & solverMatrixIdx)%solverColToEquationsSetMaps(equationsSetIdx)%SolverColToStaticEquationsMaps( &
                  !  & columnIdx)%numberOfLinearEquationsMatrices,5,5,solverMapping%solverColToEquationsColsMap( &
                  !  & solverMatrixIdx)%solverColToEquationsSetMaps(equationsSetIdx)%SolverColToStaticEquationsMaps( &
                  !  & columnIdx)%equationsMatrixNumbers,'("            Equation matrices numbers :",5(X,I13))', &
                  !  & '(36X,5(X,I13))',err,error,*999)
                  !CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,solverMapping%solverColToEquationsColsMap( &
                  !  & solverMatrixIdx)%solverColToEquationsSetMaps(equationsSetIdx)%SolverColToStaticEquationsMaps( &
                  !  & columnIdx)%numberOfLinearEquationsMatrices,5,5,solverMapping%solverColToEquationsColsMap( &
                  !  & solverMatrixIdx)%solverColToEquationsSetMaps(equationsSetIdx)%SolverColToStaticEquationsMaps( &
                  !  & columnIdx)%equationsColNumbers,'("            Equation column numbers   :",5(X,I13))', &
                  !  & '(36X,5(X,I13))',err,error,*999)
                  !CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,solverMapping%solverColToEquationsColsMap( &
                  !  & solverMatrixIdx)%solverColToEquationsSetMaps(equationsSetIdx)%SolverColToStaticEquationsMaps( &
                  !  & columnIdx)%numberOfLinearEquationsMatrices,5,5,solverMapping%solverColToEquationsColsMap( &
                  !  & solverMatrixIdx)%solverColToEquationsSetMaps(equationsSetIdx)%SolverColToStaticEquationsMaps( &
                  !  & columnIdx)%couplingCoefficients,'("            Coupling coefficients     :",5(X,E13.6))', &
                  !  & '(36X,5(X,E13.6))',err,error,*999)
                ENDIF
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"            Jacobian column number     = ", &
                  & solverMapping%solverColToEquationsColsMap(solverMatrixIdx)%solverColToEquationsSetMaps( &
                  & equationsSetIdx)%SolverColToStaticEquationsMaps(columnIdx)%jacobianColNumber,err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"            Jacobian coupling coeff    = ", &
                  & solverMapping%solverColToEquationsColsMap(solverMatrixIdx)%solverColToEquationsSetMaps( &
                  & equationsSetIdx)%SolverColToStaticEquationsMaps(columnIdx)%jacobianCouplingCoefficient,err,error,*999)
              ELSE
                 CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"            Number of static equations matrices mapped to  = ", &
                  & 0_INTG,err,error,*999)
              ENDIF
            ENDDO !columnIdx
          ENDDO !equationsSetIdx
        ENDIF
      ENDDO !solverMatrixIdx
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Solver DOF to field DOFs mappings:",err,error,*999)
      DO solverMatrixIdx=1,solverMapping%numberOfSolverMatrices
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Solver matrix : ",solverMatrixIdx,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Number of DOFs = ",solverMapping% &
          & solverColToEquationsColsMap(solverMatrixIdx)%numberOfDofs,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Total number of DOFs = ",solverMapping% &
          & solverColToEquationsColsMap(solverMatrixIdx)%totalNumberOfDofs,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Number of global DOFs = ",solverMapping% &
          & solverColToEquationsColsMap(solverMatrixIdx)%numberOfGlobalDofs,err,error,*999)
        ALLOCATE(variableTypes(solverMapping%numberOfEquationsSets+solverMapping%numberOfInterfaceConditions),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate variable types.",err,error,*999)
        DO dofIdx=1,solverMapping%solverColToEquationsColsMap(solverMatrixIdx)%totalNumberOfDofs     
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Solver local DOF : ",dofIdx,err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Number of equations mapped to     = ",solverMapping% &
            & solverColToEquationsColsMap(solverMatrixIdx)%solverDofToVariableMaps(dofIdx)% &
            & numberOfEquations,err,error,*999)
          IF(solverMapping%solverColToEquationsColsMap(solverMatrixIdx)%solverDofToVariableMaps(dofIdx)% &
            & numberOfEquations>0) THEN
            CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,solverMapping%solverColToEquationsColsMap( &
              & solverMatrixIdx)%solverDofToVariableMaps(dofIdx)%numberOfEquations,5,5,solverMapping% &
              & solverColToEquationsColsMap(solverMatrixIdx)%solverDofToVariableMaps(dofIdx)%equationsIndices, &
              & '("        Equations types       :",5(X,I13))','(31X,5(X,I13))',err,error,*999)
            CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,solverMapping%solverColToEquationsColsMap( &
              & solverMatrixIdx)%solverDofToVariableMaps(dofIdx)%numberOfEquations,5,5,solverMapping% &
              & solverColToEquationsColsMap(solverMatrixIdx)%solverDofToVariableMaps(dofIdx)%equationsIndices, &
              & '("        Equations indices     :",5(X,I13))','(31X,5(X,I13))',err,error,*999)
            DO variableIdx=1,solverMapping%solverColToEquationsColsMap(solverMatrixIdx)%solverDofToVariableMaps( &
              & dofIdx)%numberOfEquations
              variableTypes(variableIdx)=solverMapping%solverColToEquationsColsMap(solverMatrixIdx)% &
              & solverDofToVariableMaps(dofIdx)%variable(variableIdx)%ptr%VARIABLE_TYPE
            ENDDO !variableIdx
            CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,solverMapping%solverColToEquationsColsMap( &
              & solverMatrixIdx)%solverDofToVariableMaps(dofIdx)%numberOfEquations,5,5,variableTypes, &
              & '("        Variable types        :",5(X,I13))','(31X,5(X,I13))',err,error,*999)
            CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,solverMapping%solverColToEquationsColsMap( &
              & solverMatrixIdx)%solverDofToVariableMaps(dofIdx)%numberOfEquations,5,5,solverMapping% &
              & solverColToEquationsColsMap(solverMatrixIdx)%solverDofToVariableMaps(dofIdx)%variableDof, &
              & '("        Variable DOFs         :",5(X,I13))','(31X,5(X,I13))',err,error,*999)
            CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,solverMapping%solverColToEquationsColsMap( &
              & solverMatrixIdx)%solverDofToVariableMaps(dofIdx)%numberOfEquations,5,5,solverMapping% &
              & solverColToEquationsColsMap(solverMatrixIdx)%solverDofToVariableMaps(dofIdx)% &
              & variableCoefficient, & 
              & '("        Variable coefficients :",5(X,E13.6))','(31X,5(X,E13.6))',err,error,*999)
            CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,solverMapping%solverColToEquationsColsMap( &
              & solverMatrixIdx)%solverDofToVariableMaps(dofIdx)%numberOfEquations,5,5,solverMapping% &
              & solverColToEquationsColsMap(solverMatrixIdx)%solverDofToVariableMaps(dofIdx)% &
              & additiveConstant, &
              & '("        Additive constants    :",5(X,E13.6))','(31X,5(X,E13.6))',err,error,*999)
          ENDIF
        ENDDO !dofIdx
        IF(ALLOCATED(variableTypes)) DEALLOCATE(variableTypes)
      ENDDO !solverMatrixIdx
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Equation sets to solver mappings:",err,error,*999)
      DO equationsSetIdx=1,solverMapping%numberOfEquationsSets
        equationsSet=>solverMapping%equationsSets(equationsSetIdx)%ptr
        equations=>equationsSet%equations
        equationsMapping=>equations%EQUATIONS_MAPPING
        dynamicMapping=>equationsMapping%DYNAMIC_MAPPING
        linearMapping=>equationsMapping%LINEAR_MAPPING
        nonlinearMapping=>equationsMapping%NONLINEAR_MAPPING
        rhsMapping=>equationsMapping%RHS_MAPPING
        lhsMapping=>equationsMapping%lhsMapping
        sourceMapping=>equationsMapping%SOURCE_MAPPING
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Equations set index : ",equationsSetIdx,err,error,*999)
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"      Equations sets rows to solver rows mappings:",err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE, "        Number of equations set rows = ",equationsMapping% &
         & TOTAL_NUMBER_OF_ROWS,err,error,*999)
        DO rowIdx=1,equationsMapping%TOTAL_NUMBER_OF_ROWS
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Equations set row : ",rowIdx,err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"          Number of solver rows mapped to   = ",solverMapping% &
            & equationsSetToSolverMap(equationsSetIdx)%equationsRowToSolverRowsMaps(rowIdx)%numberOfSolverRows, &
            & err,error,*999)
          IF(solverMapping%equationsSetToSolverMap(equationsSetIdx)%equationsRowToSolverRowsMaps(rowIdx)% &
            & numberOfSolverRows>0) THEN
            CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,solverMapping%equationsSetToSolverMap(equationsSetIdx)% &
              & equationsRowToSolverRowsMaps(rowIdx)%numberOfSolverRows,5,5,solverMapping% &
              & equationsSetToSolverMap(equationsSetIdx)%equationsRowToSolverRowsMaps(rowIdx)%solverRows, &
              & '("          Solver row numbers    :",5(X,I13))','(33X,5(X,I13))',err,error,*999)
            CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,solverMapping%equationsSetToSolverMap(equationsSetIdx)% &
              & equationsRowToSolverRowsMaps(rowIdx)%numberOfSolverRows,5,5,solverMapping% &
              & equationsSetToSolverMap(equationsSetIdx)%equationsRowToSolverRowsMaps(rowIdx)%couplingCoefficients, &
              & '("          Coupling coefficients :",5(X,E13.6))','(33X,5(X,E13.6))',err,error,*999)
          ENDIF
        ENDDO !rowIdx
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"      Solver matrix indexing:",err,error,*999)
        DO solverMatrixIdx=1,solverMapping%numberOfSolverMatrices
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Solver matrix : ",solverMatrixIdx,err,error,*999)
          CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"          Interface conditions affecting:", &
              & err,error,*999)
           CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"            Number of interface conditions = ",solverMapping% &
            & equationsSetToSolverMap(equationsSetIdx)%numberOfInterfaceConditions,err,error,*999)
          DO interfaceConditionIdx=1,solverMapping%equationsSetToSolverMap(equationsSetIdx)%numberOfInterfaceConditions
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"            Interface condition : ",interfaceConditionIdx, &
              & err,error,*999)
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"              Interface condition index = ",solverMapping% &
              & equationsSetToSolverMap(equationsSetIdx)%equationsToSolverMatrixMapsInterface( &
              & interfaceConditionIdx)%interfaceConditionIndex,err,error,*999)
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"              Interface matrix number   = ",solverMapping% &
              & equationsSetToSolverMap(equationsSetIdx)%equationsToSolverMatrixMapsInterface( &
              & interfaceConditionIdx)%interfaceMatrixNumber,err,error,*999)
          ENDDO !interfaceConditionIdx                                 
          IF(ASSOCIATED(dynamicMapping)) THEN
           CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"          Dynamic equations matrix columns to solver matrix columns:", &
              & err,error,*999)
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"            Number of dynamic equations matrices = ",solverMapping% &
              & equationsSetToSolverMap(equationsSetIdx)%equationsToSolverMatrixMapsSm(solverMatrixIdx)% &
              & numberOfDynamicEquationsMatrices,err,error,*999)
            DO equationsMatrixIdx=1,solverMapping%equationsSetToSolverMap(equationsSetIdx)% &
              & equationsToSolverMatrixMapsSm(solverMatrixIdx)%numberOfDynamicEquationsMatrices
              equationsToSolverMap=>solverMapping%equationsSetToSolverMap(equationsSetIdx)% &
                & equationsToSolverMatrixMapsSm(solverMatrixIdx)%dynamicEquationsToSolverMatrixMaps( &
                & equationsMatrixIdx)%ptr
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"              Equations matrix index : ",equationsMatrixIdx, &
                & err,error,*999)
              equationsMatrix=equationsToSolverMap%equationsMatrixNumber
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"                Equations matrix number = ",equationsMatrix, &
                & err,error,*999)
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"                Solver matrix number    = ",equationsToSolverMap% &
                & solverMatrixNumber,err,error,*999)
              DO columnIdx=1,dynamicMapping%EQUATIONS_MATRIX_TO_VAR_MAPS(equationsMatrix)%NUMBER_OF_COLUMNS
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"                Equations matrix column : ",columnIdx, &
                  & err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"                  Number of solver columns mapped to = ", &
                  & equationsToSolverMap%equationsColToSolverColsMap(columnIdx)%numberOfSolverCols,err,error,*999)
                IF(equationsToSolverMap%equationsColToSolverColsMap(columnIdx)%numberOfSolverCols>0) THEN
                  CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,equationsToSolverMap%equationsColToSolverColsMap( &
                    & columnIdx)%numberOfSolverCols,5,5,equationsToSolverMap%equationsColToSolverColsMap(columnIdx)% &
                    & solverCols,'("                  Solver columns         :",5(X,I13))','(42X,5(X,I13))',err,error,*999)
                  CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,equationsToSolverMap%equationsColToSolverColsMap( &
                    & columnIdx)%numberOfSolverCols,5,5,equationsToSolverMap%equationsColToSolverColsMap(columnIdx)% &
                    & couplingCoefficients,'("                  Coupling coefficients  :",5(X,E13.6))','(42X,5(X,E13.6))', &
                    & err,error,*999)
                ENDIF
              ENDDO !columnIdx
            ENDDO !equationsMatrixIdx
          ELSE
            IF(ASSOCIATED(linearMapping)) THEN
                CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"          Linear equations matrix columns to solver matrix columns:", &
                & err,error,*999)
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"            Number of liner equations matrices = ",solverMapping% &
                & equationsSetToSolverMap(equationsSetIdx)%equationsToSolverMatrixMapsSm(solverMatrixIdx)% &
                & numberOfLinearEquationsMatrices,err,error,*999)
              DO equationsMatrixIdx=1,solverMapping%equationsSetToSolverMap(equationsSetIdx)% &
                & equationsToSolverMatrixMapsSm(solverMatrixIdx)%numberOfLinearEquationsMatrices
                equationsToSolverMap=>solverMapping%equationsSetToSolverMap(equationsSetIdx)% &
                  & equationsToSolverMatrixMapsSm(solverMatrixIdx)%linearEquationsToSolverMatrixMaps( &
                  & equationsMatrixIdx)%ptr
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"            Equations matrix index : ",equationsMatrixIdx, &
                  & err,error,*999)
                equationsMatrix=equationsToSolverMap%equationsMatrixNumber
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"              Equations matrix number = ",equationsMatrix, &
                  & err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"              Solver matrix number    = ", &
                  & equationsToSolverMap%solverMatrixNumber,err,error,*999)
                DO columnIdx=1,linearMapping%EQUATIONS_MATRIX_TO_VAR_MAPS(equationsMatrix)%NUMBER_OF_COLUMNS
                  CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"              Equations matrix column : ",columnIdx, &
                    & err,error,*999)
                  CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"                Number of solver columns mapped to = ", &
                    & equationsToSolverMap%equationsColToSolverColsMap(columnIdx)%numberOfSolverCols,err,error,*999)
                  IF(equationsToSolverMap%equationsColToSolverColsMap(columnIdx)%numberOfSolverCols>0) THEN
                    CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,equationsToSolverMap%equationsColToSolverColsMap( &
                      & columnIdx)%numberOfSolverCols,5,5,equationsToSolverMap%equationsColToSolverColsMap( &
                      & columnIdx)%solverCols,'("                Solver columns         :",5(X,I13))','(40X,5(X,I13))', &
                      & err,error,*999)
                    CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,equationsToSolverMap%equationsColToSolverColsMap( &
                      & columnIdx)%numberOfSolverCols,5,5,equationsToSolverMap%equationsColToSolverColsMap( &
                      & columnIdx)%couplingCoefficients, &
                      & '("                Coupling coefficients  :",5(X,E13.6))','(40X,5(X,E13.6))',err,error,*999)
                  ENDIF
                ENDDO !columnIdx
              ENDDO !equationsMatrixIdx
            ENDIF
            IF(ASSOCIATED(nonlinearMapping)) THEN
              CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"          Jacobian equations matrix columns to solver matrix columns:", &
                & err,error,*999)
              DO equationsMatrixIdx=1,nonlinearMapping%NUMBER_OF_RESIDUAL_VARIABLES
                jacobianToSolverMap=>solverMapping%equationsSetToSolverMap(equationsSetIdx)% &
                  & equationsToSolverMatrixMapsSm(solverMatrixIdx)%jacobianToSolverMatrixMaps(equationsMatrixIdx)%ptr
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"            Solver matrix number    = ",jacobianToSolverMap% &
                  & solverMatrixNumber,err,error,*999)
                DO columnIdx=1,nonlinearMapping%JACOBIAN_TO_VAR_MAP(equationsMatrixIdx)%NUMBER_OF_COLUMNS
                  CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"            Equations matrix column : ",columnIdx,err,error,*999)
                  CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"              Number of solver columns mapped to = ", &
                    & jacobianToSolverMap%jacobianColToSolverColsMap(columnIdx)%numberOfSolverCols,err,error,*999)
                  IF(jacobianToSolverMap%jacobianColToSolverColsMap(columnIdx)%numberOfSolverCols>0) THEN
                    CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,jacobianToSolverMap%jacobianColToSolverColsMap( &
                      & columnIdx)%numberOfSolverCols,5,5,jacobianToSolverMap%jacobianColToSolverColsMap(columnIdx)% &
                      & solverCols,'("              Solver columns         :",5(X,I13))','(38X,5(X,I13))',err,error,*999)
                    CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,jacobianToSolverMap%jacobianColToSolverColsMap( &
                      & columnIdx)%numberOfSolverCols,5,5,jacobianToSolverMap%jacobianColToSolverColsMap(columnIdx)% &
                      & couplingCoefficients,'("              Coupling coefficients  :",5(X,E13.6))','(38X,5(X,E13.6))', &
                      & err,error,*999)
                  ENDIF
                ENDDO !columnIdx
              ENDDO !equationsMatrixIdx
            ENDIF
          ENDIF
          CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"          Variable DOFs to solver matrix DOFs:",err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"            Number of variables = ",solverMapping% &
            & equationsSetToSolverMap(equationsSetIdx)%equationsToSolverMatrixMapsSm(solverMatrixIdx)% &
            & numberOfVariables,err,error,*999) 
          CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,solverMapping%equationsSetToSolverMap(equationsSetIdx)% &
            & equationsToSolverMatrixMapsSm(solverMatrixIdx)%numberOfVariables,5,5,solverMapping% &
            & equationsSetToSolverMap(equationsSetIdx)%equationsToSolverMatrixMapsSm(solverMatrixIdx)% &
            & variableTypes,'("            Variable types :",5(X,I13))','(28X,5(X,I13))',err,error,*999)
          DO variableIdx=1,solverMapping%equationsSetToSolverMap(equationsSetIdx)% &
            & equationsToSolverMatrixMapsSm(solverMatrixIdx)%numberOfVariables
            dependentVariable=>solverMapping%equationsSetToSolverMap(equationsSetIdx)% &
              & equationsToSolverMatrixMapsSm(solverMatrixIdx)%variables(variableIdx)%ptr
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"            Variable index : ",variableIdx,err,error,*999)
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"              Number of variable DOFs = ",dependentVariable% &
              & NUMBER_OF_DOFS,err,error,*999)
            DO localDof=1,dependentVariable%TOTAL_NUMBER_OF_DOFS
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"              Variable DOF : ",localDof,err,error,*999)
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"                Solver column number = ",solverMapping% &
                & equationsSetToSolverMap(equationsSetIdx)%equationsToSolverMatrixMapsSm(solverMatrixIdx)% &
                & variableToSolverColMaps(variableIdx)%columnNumbers(localDof),err,error,*999)
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"                Coupling coefficient = ",solverMapping% &
                & equationsSetToSolverMap(equationsSetIdx)%equationsToSolverMatrixMapsSm(solverMatrixIdx)% &
                & variableToSolverColMaps(variableIdx)%couplingCoefficients(localDof),err,error,*999)
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"                Additive constant    = ",solverMapping% &
                & equationsSetToSolverMap(equationsSetIdx)%equationsToSolverMatrixMapsSm(solverMatrixIdx)% &
                & variableToSolverColMaps(variableIdx)%additiveConstants(localDof),err,error,*999)              
            ENDDO !localDof
          ENDDO !variableIdx
        ENDDO !solverMatrixIdx
        IF(ASSOCIATED(dynamicMapping)) THEN
          CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"      Dynamic equations matrix indexing:",err,error,*999)
          DO equationsMatrixIdx=1,dynamicMapping%NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Equations matrix : ",equationsMatrixIdx,err,error,*999)
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"          Number of solver matrices = ",solverMapping% &
              & equationsSetToSolverMap(equationsSetIdx)%equationsToSolverMatrixMapsEm(equationsMatrixIdx)% &
              & numberOfSolverMatrices,err,error,*999)
            DO solverMatrixIdx=1,solverMapping%equationsSetToSolverMap(equationsSetIdx)% &
              & equationsToSolverMatrixMapsEm(equationsMatrixIdx)%numberOfSolverMatrices
              equationsToSolverMap=>solverMapping%equationsSetToSolverMap(equationsSetIdx)% &
                & equationsToSolverMatrixMapsEm(equationsMatrixIdx)%equationsToSolverMatrixMaps( &
                & solverMatrixIdx)%ptr
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"          Solver matrix index : ",solverMatrixIdx,err,error,*999)
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"            Equations matrix number = ",equationsToSolverMap% &
                & equationsMatrixNumber,err,error,*999)
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"            Solver matrix number    = ",equationsToSolverMap% &
                & solverMatrixNumber,err,error,*999)            
            ENDDO !solverMatrixIdx
          ENDDO !equationsMatrixIdx
        ELSE
          IF(ASSOCIATED(linearMapping)) THEN
            CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"      Linear equations matrix indexing:",err,error,*999)
            DO equationsMatrixIdx=1,linearMapping%NUMBER_OF_LINEAR_EQUATIONS_MATRICES
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Equations matrix : ",equationsMatrixIdx,err,error,*999)
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"          Number of solver matrices = ",solverMapping% &
                & equationsSetToSolverMap(equationsSetIdx)%equationsToSolverMatrixMapsEm(equationsMatrixIdx)% &
                & numberOfSolverMatrices,err,error,*999)
              DO solverMatrixIdx=1,solverMapping%equationsSetToSolverMap(equationsSetIdx)% &
                & equationsToSolverMatrixMapsEm(equationsMatrixIdx)%numberOfSolverMatrices
                equationsToSolverMap=>solverMapping%equationsSetToSolverMap(equationsSetIdx)% &
                  & equationsToSolverMatrixMapsEm(equationsMatrixIdx)%equationsToSolverMatrixMaps( &
                  & solverMatrixIdx)%ptr
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"          Solver matrix index : ",solverMatrixIdx,err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"            Equations matrix number = ",equationsToSolverMap% &
                  & equationsMatrixNumber,err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"            Solver matrix number    = ",equationsToSolverMap% &
                  & solverMatrixNumber,err,error,*999)            
              ENDDO !solverMatrixIdx
            ENDDO !equationsMatrixIdx
          ENDIF
          IF(ASSOCIATED(nonlinearMapping)) THEN
            CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"      Jacobian matrix indexing:",err,error,*999)
            DO equationsMatrixIdx=1,nonlinearMapping%NUMBER_OF_RESIDUAL_VARIABLES
              jacobianToSolverMap=>solverMapping%equationsSetToSolverMap(equationsSetIdx)% &
                & equationsToSolverMatrixMapsJm(equationsMatrixIdx)%ptr
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Solver matrix number    = ",jacobianToSolverMap% &
                & solverMatrixNumber,err,error,*999)
            ENDDO
          ENDIF
        ENDIF
      ENDDO !equationsSetIdx
      IF(solverMapping%numberOfInterfaceConditions>0) THEN
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Interface conditions to solver mappings:",err,error,*999)
        DO interfaceConditionIdx=1,solverMapping%numberOfInterfaceConditions
          interfaceCondition=>solverMapping%interfaceConditions(interfaceConditionIdx)%ptr
          interfaceEquations=>interfaceCondition%INTERFACE_EQUATIONS
          interfaceMapping=>interfaceEquations%INTERFACE_MAPPING
          CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"    Interface to equations sets mapping:",err,error,*999)
          DO solverMatrixIdx=1,solverMapping%numberOfSolverMatrices
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Solver matrix : ",solverMatrixIdx,err,error,*999)
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Number of equations sets = ",solverMapping% &
              & interfaceConditionToSolverMap(interfaceConditionIdx)%numberOfEquationsSets,err,error,*999)
            DO equationsSetIdx=1,solverMapping%interfaceConditionToSolverMap(interfaceConditionIdx)% &
              & numberOfEquationsSets
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Equations set : ",equationsSetIdx,err,error,*999)
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"          Equations set index     = ",solverMapping% &
                & interfaceConditionToSolverMap(interfaceConditionIdx)%interfaceToSolverMatrixMapsEquations( &
                & equationsSetIdx)%equationsSetIndex,err,error,*999)
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"          Interface matrix number = ",solverMapping% &
                & interfaceConditionToSolverMap(interfaceConditionIdx)%interfaceToSolverMatrixMapsEquations( &
                & equationsSetIdx)%interfaceMatrixIndex,err,error,*999)
            ENDDO !equationsSetIdx
          ENDDO !solverMatrixIdx
          CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"    Interface condition rows to solver rows mappings:",err,error,*999)
          DO interfaceMatrixIdx=1,interfaceMapping%NUMBER_OF_INTERFACE_MATRICES
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Interface matrix idx : ",interfaceMatrixIdx,err,error,*999)
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Number of interface rows = ",interfaceMapping% &
              & INTERFACE_MATRIX_ROWS_TO_VAR_MAPS(interfaceMatrixIdx)%TOTAL_NUMBER_OF_ROWS,err,error,*999)
            DO rowIdx=1,interfaceMapping%INTERFACE_MATRIX_ROWS_TO_VAR_MAPS(interfaceMatrixIdx)%TOTAL_NUMBER_OF_ROWS
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Interface row : ",rowIdx,err,error,*999)
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"          Number of solver rows mapped to = ",solverMapping% &
                & interfaceConditionToSolverMap(interfaceConditionIdx)%interfaceToSolverMatrixMapsIm( &
                & interfaceMatrixIdx)%interfaceRowToSolverRowsMap(rowIdx)%numberOfSolverRows,err,error,*999)
              IF(solverMapping%interfaceConditionToSolverMap(interfaceConditionIdx)%interfaceToSolverMatrixMapsIm( &
                & interfaceMatrixIdx)%interfaceRowToSolverRowsMap(rowIdx)%numberOfSolverRows>0) THEN
                CALL WRITE_STRING_FMT_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"          Solver row numbers    : ",solverMapping% &
                  & interfaceConditionToSolverMap(interfaceConditionIdx)%interfaceToSolverMatrixMapsIm( &
                  & interfaceMatrixIdx)%interfaceRowToSolverRowsMap(rowIdx)%solverRow,"(I13)",err,error,*999)
                CALL WRITE_STRING_FMT_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"          Coupling coefficients : ",solverMapping% &
                  & interfaceConditionToSolverMap(interfaceConditionIdx)%interfaceToSolverMatrixMapsIm( &
                  & interfaceMatrixIdx)%interfaceRowToSolverRowsMap(rowIdx)%couplingCoefficient,"(E13.6)",err,error,*999)
              ENDIF
            ENDDO !rowIdx
          ENDDO !interfaceMatrixIdx
          CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"    Interface condition column to solver rows mappings:", &
            & err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Number of interface condition columns = ",interfaceMapping% &
            & TOTAL_NUMBER_OF_COLUMNS,err,error,*999)
          DO columnIdx=1,interfaceMapping%TOTAL_NUMBER_OF_COLUMNS
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Interface condition column : ",columnIdx,err,error,*999)
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Number of rows mapped to = ",solverMapping% &
              & interfaceConditionToSolverMap(interfaceConditionIdx)%interfaceColumnToSolverRowsMaps(columnIdx)% &
              & numberOfSolverRows,err,error,*999)
            IF(solverMapping%interfaceConditionToSolverMap(interfaceConditionIdx)% &
              & interfaceColumnToSolverRowsMaps(columnIdx)%numberOfSolverRows>0) THEN
              CALL WRITE_STRING_FMT_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Solver row number    : ",solverMapping% &
                & interfaceConditionToSolverMap(interfaceConditionIdx)%interfaceColumnToSolverRowsMaps(columnIdx)% &
                & solverRow, "(I13)",err,error,*999)
              CALL WRITE_STRING_FMT_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Coupling coefficient : ",solverMapping% &
                & interfaceConditionToSolverMap(interfaceConditionIdx)%interfaceColumnToSolverRowsMaps(columnIdx)% &
                & couplingCoefficient, "(E13.6)",err,error,*999)
            ENDIF
          ENDDO !columnIdx
          CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"    Solver matrix indexing:",err,error,*999)
          DO solverMatrixIdx=1,solverMapping%numberOfSolverMatrices
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Solver matrix : ",solverMatrixIdx,err,error,*999)        
            CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"        Interface equations matrix rows to solver matrix columns:", &
              & err,error,*999)
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"          Number of interface matrices = ",solverMapping% &
              & interfaceConditionToSolverMap(interfaceConditionIdx)%interfaceToSolverMatrixMapsSm(solverMatrixIdx)% &
              & numberOfInterfaceMatrices,err,error,*999)
            DO interfaceMatrixIdx=1,solverMapping%interfaceConditionToSolverMap(interfaceConditionIdx)% &
              & interfaceToSolverMatrixMapsSm(solverMatrixIdx)%numberOfInterfaceMatrices
              interfaceToSolverMap=>solverMapping%interfaceConditionToSolverMap(interfaceConditionIdx)% &
                & interfaceToSolverMatrixMapsSm(solverMatrixIdx)%interfaceEquationsToSolverMatrixMaps( &
                & interfaceMatrixIdx)%ptr
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"          Interface matrix index : ",interfaceMatrixIdx, &
                & err,error,*999)
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"            Interface matrix number = ",interfaceToSolverMap% &
                & interfaceMatrixNumber,err,error,*999)
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"            Solver matrix number    = ",interfaceToSolverMap% &
                & solverMatrixNumber,err,error,*999)
              DO rowIdx=1,interfaceMapping%INTERFACE_MATRIX_ROWS_TO_VAR_MAPS(interfaceMatrixIdx)%NUMBER_OF_ROWS
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"            Interface matrix row : ",rowIdx, &
                  & err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"              Number of solver columns mapped to = ", &
                  & interfaceToSolverMap%interfaceRowToSolverColsMap(rowIdx)%numberOfSolverCols,err,error,*999)
                IF(interfaceToSolverMap%interfaceRowToSolverColsMap(rowIdx)%numberOfSolverCols>0) THEN
                  CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,interfaceToSolverMap%interfaceRowToSolverColsMap( &
                    & rowIdx)%numberOfSolverCols,5,5,interfaceToSolverMap%interfaceRowToSolverColsMap(rowIdx)% &
                    & solverCols,'("              Solver columns         :",5(X,I13))','(38X,5(X,I13))',err,error,*999)
                  CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,interfaceToSolverMap%interfaceRowToSolverColsMap( &
                    & rowIdx)%numberOfSolverCols,5,5,interfaceToSolverMap%interfaceRowToSolverColsMap(rowIdx)% &
                    & couplingCoefficients,'("              Coupling coefficients  :",5(X,E13.6))','(38X,5(X,E13.6))', &
                    & err,error,*999)
                ENDIF
              ENDDO !rowIdx
            ENDDO !interfaceMatrixIdx
            CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"        Variable dofs to solver matrix dofs:",err,error,*999)
            CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"          Lagrange variables:",err,error,*999)
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"            Lagrange variable type = ",solverMapping% &
              & interfaceConditionToSolverMap(interfaceConditionIdx)%interfaceToSolverMatrixMapsSm(solverMatrixIdx)% &
              & lagrangeVariableType,err,error,*999)
            lagrangeVariable=>solverMapping%interfaceConditionToSolverMap(interfaceConditionIdx)% &
              & interfaceToSolverMatrixMapsSm(solverMatrixIdx)%lagrangeVariable
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"            Number of Lagrange variable dofs = ",lagrangeVariable% &
              & NUMBER_OF_DOFS,err,error,*999)
            DO localDof=1,lagrangeVariable%TOTAL_NUMBER_OF_DOFS
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"            Variable dof : ",localDof,err,error,*999)
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"              Solver column number = ",solverMapping% &
                & interfaceConditionToSolverMap(interfaceConditionIdx)%interfaceToSolverMatrixMapsSm( &
                & solverMatrixIdx)%lagrangeVariableToSolverColMap%columnNumbers(localDof),err,error,*999)
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"              Coupling coefficient = ",solverMapping% &
                & interfaceConditionToSolverMap(interfaceConditionIdx)%interfaceToSolverMatrixMapsSm( &
                & solverMatrixIdx)%lagrangeVariableToSolverColMap%couplingCoefficients(localDof),err,error,*999)
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"              Additive constant    = ",solverMapping% &
                & interfaceConditionToSolverMap(interfaceConditionIdx)%interfaceToSolverMatrixMapsSm( &
                & solverMatrixIdx)%lagrangeVariableToSolverColMap%additiveConstants(localDof),err,error,*999)              
            ENDDO !localDof
            CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"          Dependent variables:",err,error,*999)
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"            Number of dependent variables = ",solverMapping% &
              & interfaceConditionToSolverMap(interfaceConditionIdx)%interfaceToSolverMatrixMapsSm(solverMatrixIdx)% &
              & numberOfDependentVariables,err,error,*999)
            CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,solverMapping%interfaceConditionToSolverMap( &
              & interfaceConditionIdx)%interfaceToSolverMatrixMapsSm(solverMatrixIdx)%numberOfDependentVariables, &
              & 5,5,solverMapping%interfaceConditionToSolverMap(interfaceConditionIdx)%interfaceToSolverMatrixMapsSm( &
              & solverMatrixIdx)%dependentVariableTypes,'("            Dependent variable types :",5(X,I13))', &
              & '(38X,5(X,I13))',err,error,*999) 
            DO variableIdx=1,solverMapping%interfaceConditionToSolverMap(interfaceConditionIdx)% &
              & interfaceToSolverMatrixMapsSm(solverMatrixIdx)%numberOfDependentVariables
              dependentVariable=>solverMapping%interfaceConditionToSolverMap(interfaceConditionIdx)% &
                & interfaceToSolverMatrixMapsSm(solverMatrixIdx)%dependentVariables(variableIdx)%ptr
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"            Dependent variable index : ",variableIdx,err,error,*999)
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"              Number of dependent variable dofs = ", &
                & dependentVariable%NUMBER_OF_DOFS,err,error,*999)
              DO localDof=1,dependentVariable%TOTAL_NUMBER_OF_DOFS
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"              Variable dof : ",localDof,err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"                Solver column number = ",solverMapping% &
                  & interfaceConditionToSolverMap(interfaceConditionIdx)%interfaceToSolverMatrixMapsSm( &
                  & solverMatrixIdx)%dependentVariableToSolverColMaps(variableIdx)%columnNumbers(localDof), &
                  & err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"                Coupling coefficient = ",solverMapping% &
                  & interfaceConditionToSolverMap(interfaceConditionIdx)%interfaceToSolverMatrixMapsSm( &
                  & solverMatrixIdx)%dependentVariableToSolverColMaps(variableIdx)%couplingCoefficients(localDof), &
                  & err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"                Additive constant    = ",solverMapping% &
                  & interfaceConditionToSolverMap(interfaceConditionIdx)%interfaceToSolverMatrixMapsSm( &
                  & solverMatrixIdx)%dependentVariableToSolverColMaps(variableIdx)%additiveConstants(localDof), &
                  & err,error,*999)              
              ENDDO !localDof
            ENDDO !variableIdx
          ENDDO !solverMatrixIdx        
          CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"    Interface equations matrix indexing:",err,error,*999)
          DO interfaceMatrixIdx=1,interfaceMapping%NUMBER_OF_INTERFACE_MATRICES
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Interface matrix : ",interfaceMatrixIdx,err,error,*999)
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Number of solver matrices = ",solverMapping% &
              & interfaceConditionToSolverMap(interfaceConditionIdx)%interfaceToSolverMatrixMapsIm( &
              & interfaceMatrixIdx)%numberOfSolverMatrices,err,error,*999)
            DO solverMatrixIdx=1,solverMapping%interfaceConditionToSolverMap(interfaceConditionIdx)% &
              & interfaceToSolverMatrixMapsIm(interfaceMatrixIdx)%numberOfSolverMatrices
              interfaceToSolverMap=>solverMapping%interfaceConditionToSolverMap(interfaceConditionIdx)% &
                & interfaceToSolverMatrixMapsIm(interfaceMatrixIdx)%interfaceToSolverMatrixMaps( &
                & solverMatrixIdx)%ptr
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Solver matrix index : ",solverMatrixIdx,err,error,*999)
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"          Interface matrix number = ",interfaceToSolverMap% &
                & interfaceMatrixNumber,err,error,*999)
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"          Solver matrix number    = ",interfaceToSolverMap% &
                & solverMatrixNumber,err,error,*999)            
            ENDDO !solverMatrixIdx
          ENDDO !equationsMatrixIdx
          CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"    Interface column to solver rows mapping:",err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Number of columns = ",interfaceMapping%NUMBER_OF_COLUMNS, &
            & err,error,*999)            
          DO columnIdx=1,interfaceMapping%NUMBER_OF_COLUMNS
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Column : ",columnIdx, err,error,*999)
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Number of solver rows = ",solverMapping% &
              & interfaceConditionToSolverMap(interfaceConditionIdx)%interfaceColumnToSolverRowsMaps(columnIdx)% &
              & numberOfSolverRows,err,error,*999)
            IF(solverMapping%interfaceConditionToSolverMap(interfaceConditionIdx)%interfaceColumnToSolverRowsMaps( &
              & columnIdx)%numberOfSolverRows>0) THEN
             CALL WriteStringFmtValue(DIAGNOSTIC_OUTPUT_TYPE,"        Solver row             : ",solverMapping% &
                & interfaceConditionToSolverMap(interfaceConditionIdx)%interfaceColumnToSolverRowsMaps( &
                & columnIdx)%solverRow,"(I13)",err,error,*999)
             CALL WriteStringFmtValue(DIAGNOSTIC_OUTPUT_TYPE,"        Coupling coefficients  : ",solverMapping% &
                & interfaceConditionToSolverMap(interfaceConditionIdx)%interfaceColumnToSolverRowsMaps( &
                & columnIdx)%couplingCoefficient,"(E13.6)", &
                & err,error,*999)
            ENDIF
          ENDDO !columnIdx
        ENDDO !interfaceConditionIdx
      ENDIF
    ENDIF
    
    CALL Exits("SolverMapping_Calculate")
    RETURN
999 IF(ALLOCATED(subMatrixInformation)) DEALLOCATE(subMatrixInformation)
    IF(ALLOCATED(subMatrixList)) DEALLOCATE(subMatrixList)
    IF(ALLOCATED(variablesList)) DEALLOCATE(variablesList)
    IF(ALLOCATED(variableProcessed)) DEALLOCATE(variableProcessed)
    IF(ALLOCATED(variableRankProcessed)) DEALLOCATE(variableRankProcessed)
    IF(ALLOCATED(numberOfVariableGlobalSolverDofs)) DEALLOCATE(numberOfVariableGlobalSolverDofs)
    IF(ALLOCATED(numberOfVariableLocalSolverDofs)) DEALLOCATE(numberOfVariableLocalSolverDofs)
    IF(ALLOCATED(totalNumberOfVariableLocalSolverDofs)) DEALLOCATE(totalNumberOfVariableLocalSolverDofs)    
    IF(ALLOCATED(dummyDofCoupling%globalDofs)) DEALLOCATE(dummyDofCoupling%globalDofs)
    IF(ALLOCATED(dummyDofCoupling%localDofs)) DEALLOCATE(dummyDofCoupling%localDofs)
    IF(ALLOCATED(dummyDofCoupling%coefficients)) DEALLOCATE(dummyDofCoupling%coefficients)
    CALL SolverDofCouplings_Finalise(rowCouplings,dummyErr,dummyError,*998)
998 CALL SolverDofCouplings_Finalise(columnCouplings,dummyErr,dummyError,*997)
997 CALL Errors("SolverMapping_Calculate",err,error)
    CALL Exits("SolverMapping_Calculate")
    RETURN 1
  END SUBROUTINE SolverMapping_Calculate

  !
  !================================================================================================================================
  !

  !>Finishes the process of creating a solver mapping
  SUBROUTINE SolverMapping_CreateFinish(solverMapping,err,error,*)

    !Argument variables
    TYPE(SolverMappingType), POINTER :: solverMapping !<A pointer to the solver mapping
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    CALL Enters("SolverMapping_CreateFinish",err,error,*998)

    IF(ASSOCIATED(solverMapping)) THEN
      IF(solverMapping%solverMappingFinished) THEN
        CALL FlagError("Solver mapping has already been finished",err,error,*998)
      ELSE
        IF(ASSOCIATED(solverMapping%createValuesCache)) THEN
          CALL SolverMapping_Calculate(solverMapping,err,error,*999)
          CALL SolverMapping_CreateValuesCacheFinalise(solverMapping%createValuesCache,err,error,*999)
          solverMapping%solverMappingFinished=.TRUE.            
        ELSE
          CALL FlagError("Solver mapping create values cache is not associated",err,error,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Solver mapping is not associated",err,error,*998)
    ENDIF
       
    CALL Exits("SolverMapping_CreateFinish")
    RETURN
999 CALL SolverMapping_Finalise(solverMapping,dummyErr,dummyError,*998)
998 CALL Errors("SolverMapping_CreateFinish",err,error)
    CALL Exits("SolverMapping_CreateFinish")
    RETURN 1
  END SUBROUTINE SolverMapping_CreateFinish

  !
  !================================================================================================================================
  !

  !>Starts the process of creating a solver mapping for a problem solver
  SUBROUTINE SolverMapping_CreateStart(solverEquations,solverMapping,err,error,*)

    !Argument variables
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: solverEquations !<A pointer to the solver equations to create the solver mapping on.
    TYPE(SolverMappingType), POINTER :: solverMapping !<On return, a pointer to the solver mapping. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    CALL Enters("SolverMapping_CreateStart",err,error,*999)

    IF(ASSOCIATED(solverEquations)) THEN
      IF(solverEquations%SOLVER_EQUATIONS_FINISHED) THEN
        CALL FlagError("Solver equations has already been finished.",err,error,*999)
      ELSE
        IF(ASSOCIATED(solverMapping)) THEN
          CALL FlagError("Solver mapping is already associated.",err,error,*999)
        ELSE
          NULLIFY(solverMapping)
          CALL SolverMapping_Initialise(solverEquations,err,error,*999)
          solverMapping=>solverEquations%SOLVER_MAPPING
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Solver is not associated",err,error,*999)
    ENDIF
       
    CALL Exits("SolverMapping_CreateStart")
    RETURN
999 CALL Errors("SolverMapping_CreateStart",err,error)
    CALL Exits("SolverMapping_CreateStart")
    RETURN 1
  END SUBROUTINE SolverMapping_CreateStart

  !
  !================================================================================================================================
  !

  !>Finalises a solver mapping create values cache and deallocates all memory
  SUBROUTINE SolverMapping_CreateValuesCacheFinalise(createValuesCache,err,error,*)

    !Argument variables
    TYPE(SolverMappingCreateValuesCacheType), POINTER :: createValuesCache !<A pointer to the create values cache to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: equationsSetIdx,solverMatrixIdx

    CALL Enters("SolverMapping_CreateValuesCacheFinalise",err,error,*999)

    IF(ASSOCIATED(createValuesCache)) THEN
      IF(ASSOCIATED(createValuesCache%equationsRowVariablesList)) &
        & CALL List_Destroy(createValuesCache%equationsRowVariablesList,err,error,*999)
      IF(ASSOCIATED(createValuesCache%equationsColVariablesList)) THEN
        DO solverMatrixIdx=1,SIZE(createValuesCache%equationsColVariablesList,1)
          IF(ASSOCIATED(createValuesCache%equationsColVariablesList(solverMatrixIdx)%PTR)) & 
            & CALL List_Destroy(createValuesCache%equationsColVariablesList(solverMatrixIdx)%PTR,err,error,*999)
        ENDDO !solverMatrixIdx
        DEALLOCATE(createValuesCache%equationsColVariablesList)
      ENDIF
      IF(ASSOCIATED(createValuesCache%interfaceRowVariablesList)) &
        & CALL List_Destroy(createValuesCache%interfaceRowVariablesList,err,error,*999)
      IF(ASSOCIATED(createValuesCache%interfaceColVariablesList)) THEN
        DO solverMatrixIdx=1,SIZE(createValuesCache%interfaceColVariablesList,1)
          IF(ASSOCIATED(createValuesCache%interfaceColVariablesList(solverMatrixIdx)%ptr)) & 
            & CALL List_Destroy(createValuesCache%interfaceColVariablesList(solverMatrixIdx)%PTR,err,error,*999)
        ENDDO !solverMatrixIdx
        DEALLOCATE(createValuesCache%interfaceColVariablesList)
      ENDIF
      IF(ASSOCIATED(createValuesCache%interfaceIndices)) THEN
        DO equationsSetIdx=1,SIZE(createValuesCache%interfaceIndices,1)
          IF(ASSOCIATED(createValuesCache%interfaceIndices(equationsSetIdx)%PTR)) &
            & CALL List_Destroy(createValuesCache%interfaceIndices(equationsSetIdx)%PTR, &
            & err,error,*999)
        ENDDO !equationsSetIdx
        DEALLOCATE(createValuesCache%interfaceIndices)
      ENDIF
      IF(ALLOCATED(createValuesCache%dynamicVariableType)) DEALLOCATE(createValuesCache%dynamicVariableType)
      IF(ALLOCATED(createValuesCache%matrixVariableTypes)) DEALLOCATE(createValuesCache%matrixVariableTypes)
      IF(ALLOCATED(createValuesCache%residualVariableTypes)) DEALLOCATE(createValuesCache%residualVariableTypes)
      IF(ALLOCATED(createValuesCache%lhsVariableType)) DEALLOCATE(createValuesCache%lhsVariableType)
      IF(ALLOCATED(createValuesCache%rhsVariableType)) DEALLOCATE(createValuesCache%rhsVariableType)
      IF(ALLOCATED(createValuesCache%sourceVariableType)) DEALLOCATE(createValuesCache%sourceVariableType)
      DEALLOCATE(createValuesCache)
    ENDIF
       
    CALL Exits("SolverMapping_CreateValuesCacheFinalise")
    RETURN
999 CALL Errors("SolverMapping_CreateValuesCacheFinalise",err,error)
    CALL Exits("SolverMapping_CreateValuesCacheFinalise")
    RETURN 1
  END SUBROUTINE SolverMapping_CreateValuesCacheFinalise

  !
  !================================================================================================================================
  !

  !>Initialises a solver mapping create values cache 
  SUBROUTINE SolverMapping_CreateValuesCacheInitialise(solverMapping,err,error,*)

    !Argument variables
    TYPE(SolverMappingType), POINTER :: solverMapping !<A pointer to the create values cache
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr,equationsSetIdx,solverMatrixIdx
    TYPE(VARYING_STRING) :: dummyError

    CALL Enters("SolverMapping_CreateValuesCacheInitialise",err,error,*998)

    IF(ASSOCIATED(solverMapping)) THEN
      IF(ASSOCIATED(solverMapping%createValuesCache)) THEN
        CALL FlagError("Solver mapping create values cache is already associated.",err,error,*998)
      ELSE
        ALLOCATE(solverMapping%createValuesCache,STAT=err)
        NULLIFY(solverMapping%createValuesCache%equationsRowVariablesList)
        IF(err/=0) CALL FlagError("Could not allocate solver mapping create values cache.",err,error,*999)
        ALLOCATE(solverMapping%createValuesCache%equationsColVariablesList(solverMapping%numberOfSolverMatrices),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate solver mapping create values cache equations variable list.",err,error,*999)
        NULLIFY(solverMapping%createValuesCache%interfaceRowVariablesList)
        ALLOCATE(solverMapping%createValuesCache%interfaceColVariablesList(solverMapping%numberOfSolverMatrices),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate solver mapping create values cache interface variable list.",err,error,*999)
        ALLOCATE(solverMapping%createValuesCache%interfaceIndices(solverMapping%numberOfEquationsSets),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate solver mapping create values cache interface condition indices.", &
          & err,error,*999)
        ALLOCATE(solverMapping%createValuesCache%dynamicVariableType(solverMapping%numberOfEquationsSets),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate solver mapping create values cache dynamic variable type.",err,error,*999)
        ALLOCATE(solverMapping%createValuesCache%matrixVariableTypes(0:FIELD_NUMBER_OF_VARIABLE_TYPES, &
          & solverMapping%numberOfEquationsSets,solverMapping%numberOfSolverMatrices),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate solver mapping create values cache matrix variable types.", err,error,*999)
        ALLOCATE(solverMapping%createValuesCache%residualVariableTypes(0:FIELD_NUMBER_OF_VARIABLE_TYPES, &
          & solverMapping%numberOfEquationsSets),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate solver mapping create values cache residual variable type.",err,error,*999)
        ALLOCATE(solverMapping%createValuesCache%lhsVariableType(solverMapping%numberOfEquationsSets),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate solver mapping create values cache LHS variable type.",err,error,*999)
        ALLOCATE(solverMapping%createValuesCache%rhsVariableType(solverMapping%numberOfEquationsSets),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate solver mapping create values cache RHS variable type.",err,error,*999)
        ALLOCATE(solverMapping%createValuesCache%sourceVariableType(solverMapping%numberOfEquationsSets),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate solver mapping create values cache source variable type.",err,error,*999)
        CALL List_CreateStart(solverMapping%createValuesCache%equationsRowVariablesList,err,error,*999)
        CALL List_DataTypeSet(solverMapping%createValuesCache%equationsRowVariablesList,LIST_INTG_TYPE,err,error,*999)
        CALL List_DataDimensionSet(solverMapping%createValuesCache%equationsRowVariablesList,2,err,error,*999)
        CALL List_CreateFinish(solverMapping%createValuesCache%equationsRowVariablesList,err,error,*999)
        CALL List_CreateStart(solverMapping%createValuesCache%interfaceRowVariablesList,err,error,*999)
        CALL List_DataTypeSet(solverMapping%createValuesCache%interfaceRowVariablesList,LIST_INTG_TYPE,err,error,*999)
        CALL List_DataDimensionSet(solverMapping%createValuesCache%interfaceRowVariablesList,2,err,error,*999)
        CALL List_CreateFinish(solverMapping%createValuesCache%interfaceRowVariablesList,err,error,*999)
        DO solverMatrixIdx=1,solverMapping%numberOfSolverMatrices
          NULLIFY(solverMapping%createValuesCache%equationsColVariablesList(solverMatrixIdx)%ptr)
          CALL List_CreateStart(solverMapping%createValuesCache%equationsColVariablesList(solverMatrixIdx)%ptr,err,error,*999)
          CALL List_DataTypeSet(solverMapping%createValuesCache%equationsColVariablesList(solverMatrixIdx)%ptr,LIST_INTG_TYPE, &
            & err,error,*999)
          CALL List_DataDimensionSet(solverMapping%createValuesCache%equationsColVariablesList(solverMatrixIdx)%ptr,2, &
            & err,error,*999)
          CALL List_CreateFinish(solverMapping%createValuesCache%equationsColVariablesList(solverMatrixIdx)%ptr,err,error,*999)
          NULLIFY(solverMapping%createValuesCache%interfaceColVariablesList(solverMatrixIdx)%ptr)
          CALL List_CreateStart(solverMapping%createValuesCache%interfaceColVariablesList(solverMatrixIdx)%ptr,err,error,*999)
          CALL List_DataTypeSet(solverMapping%createValuesCache%interfaceColVariablesList(solverMatrixIdx)%ptr,LIST_INTG_TYPE, &
            & err,error,*999)
          CALL List_DataDimensionSet(solverMapping%createValuesCache%interfaceColVariablesList(solverMatrixIdx)%ptr,2, &
            & err,error,*999)
          CALL List_CreateFinish(solverMapping%createValuesCache%interfaceColVariablesList(solverMatrixIdx)%ptr,err,error,*999)
        ENDDO !solver_idx
        DO equationsSetIdx=1,solverMapping%numberOfEquationsSets
          NULLIFY(solverMapping%createValuesCache%interfaceIndices(equationsSetIdx)%ptr)
          CALL List_CreateStart(solverMapping%createValuesCache%interfaceIndices(equationsSetIdx)%ptr,err,error,*999)
          CALL List_DataTypeSet(solverMapping%createValuesCache%interfaceIndices(equationsSetIdx)%ptr,LIST_INTG_TYPE,err,error,*999)
          CALL List_DataDimensionSet(solverMapping%createValuesCache%interfaceIndices(equationsSetIdx)%ptr,2,err,error,*999)
          CALL List_KeyDimensionSet(solverMapping%createValuesCache%interfaceIndices(equationsSetIdx)%ptr,1,err,error,*999)
          CALL List_CreateFinish(solverMapping%createValuesCache%interfaceIndices(equationsSetIdx)%ptr,err,error,*999)            
        ENDDO !equationsSetIdx
        solverMapping%createValuesCache%dynamicVariableType=0
        solverMapping%createValuesCache%matrixVariableTypes=0
        solverMapping%createValuesCache%residualVariableTypes=0
        solverMapping%createValuesCache%lhsVariableType=0
        solverMapping%createValuesCache%rhsVariableType=0
        solverMapping%createValuesCache%sourceVariableType=0
      ENDIF
    ELSE
      CALL FlagError("Solver mapping is not associated.",err,error,*998)
    ENDIF
       
    CALL Exits("SolverMapping_CreateValuesCacheInitialise")
    RETURN
999 CALL SolverMapping_CreateValuesCacheFinalise(solverMapping%createValuesCache,dummyErr,dummyError,*998)
998 CALL Errors("SolverMapping_CreateValuesCacheInitialise",err,error)
    CALL Exits("SolverMapping_CreateValuesCacheInitialise")
    RETURN 1
  END SUBROUTINE SolverMapping_CreateValuesCacheInitialise

  !
  !================================================================================================================================
  !

  !>Adds a variable type from an equations set dependent field to the list of variables for a particular equations variable list of a solver mapping.
  SUBROUTINE SolverMapping_CreateValuesCacheEquationVarListAdd(solverMapping,equationVarList,equationsSetIdx,variableType, &
    & err,error,*)

    !Argument variables
    TYPE(SolverMappingType), POINTER :: solverMapping !<A pointer the solver mapping to add to the var list for.
    TYPE(LIST_TYPE), POINTER :: equationVarList !<The equation variable list to add to
    INTEGER(INTG), INTENT(IN) :: equationsSetIdx !<The equations set index of the variable to add
    INTEGER(INTG), INTENT(IN) :: variableType !<The variable type to add
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: equationsSetIdx2,numberOfVariables,variableIdx,variableItem(2)
    LOGICAL :: variableFound
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet,varEquationsSet
    TYPE(FIELD_TYPE), POINTER :: dependentField,varDependentField
    TYPE(VARYING_STRING) :: localError

    CALL Enters("SolverMapping_CreateValuesCacheEquationVarListAdd",err,error,*999)

    IF(ASSOCIATED(solverMapping)) THEN
      IF(ASSOCIATED(equationVarList)) THEN
        IF(equationsSetIdx>0.AND.equationsSetIdx<=solverMapping%numberOfEquationsSets) THEN
          equationsSet=>solverMapping%equationsSets(equationsSetIdx)%ptr
          IF(ASSOCIATED(equationsSet)) THEN
            dependentField=>equationsSet%dependent%DEPENDENT_FIELD
            IF(ASSOCIATED(dependentField)) THEN
              IF(variableType/=0) THEN
                variableFound=.FALSE.
                CALL List_NumberOfItemsGet(equationVarList,numberOfVariables,err,error,*999)
                DO variableIdx=1,numberOfVariables
                  CALL List_ItemGet(equationVarList,variableIdx,variableItem,err,error,*999)
                  equationsSetIdx2=variableItem(1)
                  varEquationsSet=>solverMapping%equationsSets(equationsSetIdx2)%ptr
                  IF(ASSOCIATED(varEquationsSet)) THEN
                    varDependentField=>varEquationsSet%dependent%DEPENDENT_FIELD
                    IF(ASSOCIATED(varDependentField)) THEN
                      IF(ASSOCIATED(dependentField,varDependentField)) THEN
                        IF(variableType==variableItem(2)) THEN
                          variableFound=.TRUE.
                          EXIT
                        ENDIF
                      ENDIF
                    ELSE
                      CALL FlagError("Variable dependent field is not associated.",err,error,*999)
                    ENDIF
                  ELSE
                    CALL FlagError("Variable equations set is not associated.",err,error,*999)
                  ENDIF
                ENDDO !variableIdx
                IF(.NOT.variableFound) THEN
                  variableItem(1)=equationsSetIdx
                  variableItem(2)=variableType
                  CALL List_ItemAdd(equationVarList,variableItem,err,error,*999)
                ENDIF
              ENDIF
            ELSE
              CALL FlagError("Dependent field is not associated.",err,error,*999)
            ENDIF
          ELSE
            CALL FlagError("Equations set is not associated.",err,error,*999)
          ENDIF
        ELSE
          localError="The specified equations set index of "//TRIM(NumberToVString(equationsSetIdx,"*",err,error))// &
            & " is invalid. The index must be > 0 and <= "// &
            & TRIM(NumberToVString(solverMapping%numberOfEquationsSets,"*",err,error))//"."        
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Equation variable list is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Solver mapping is not associated.",err,error,*999)
    ENDIF
    
    CALL Exits("SolverMapping_CreateValuesCacheEquationVarListAdd")
    RETURN
999 CALL Errors("SolverMapping_CreateValuesCacheEquationVarListAdd",err,error)    
    CALL Exits("SolverMapping_CreateValuesCacheEquationVarListAdd")
    RETURN 1
   
  END SUBROUTINE SolverMapping_CreateValuesCacheEquationVarListAdd

  !
  !================================================================================================================================
  !

  !>Adds a variable type from an interface condition Lagrange field to the list of variables for a particular interface variable list of a solver mapping.
  SUBROUTINE SolverMapping_CreateValuesCacheInterfaceVarListAdd(solverMapping,interfaceVarList,interfaceConditionIdx, &
    & variableType,err,error,*)
    
    !Argument variables
    TYPE(SolverMappingType), POINTER :: solverMapping !<A pointer the solver mapping to add to the var list for.
    TYPE(LIST_TYPE), POINTER :: interfaceVarList !<The interface variables list to add to
    INTEGER(INTG), INTENT(IN) :: interfaceConditionIdx !<The interface condition index of the variable to add
    INTEGER(INTG), INTENT(IN) :: variableType !<The variable type to add
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: interfaceConditionIdx2,numberOfVariables,variableIdx,variableItem(2)
    LOGICAL :: variableFound
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: interfaceCondition,varInterfaceCondition
    TYPE(FIELD_TYPE), POINTER :: lagrangeField,varLagrangeField
    TYPE(VARYING_STRING) :: localError

    CALL Enters("SolverMapping_CreateValuesCacheInterfaceVarListAdd",err,error,*999)

    IF(ASSOCIATED(solverMapping)) THEN
      IF(ASSOCIATED(interfaceVarList)) THEN
        IF(interfaceConditionIdx>0.AND.interfaceConditionIdx<=solverMapping%numberOfInterfaceConditions) THEN
          interfaceCondition=>solverMapping%interfaceConditions(interfaceConditionIdx)%ptr
          IF(ASSOCIATED(interfaceCondition)) THEN
            SELECT CASE(interfaceCondition%method)
            CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD)
              IF(ASSOCIATED(interfaceCondition%lagrange)) THEN
                lagrangeField=>interfaceCondition%lagrange%LAGRANGE_FIELD
                IF(ASSOCIATED(lagrangeField)) THEN
                  IF(variableType/=0) THEN
                    variableFound=.FALSE.
                    CALL List_NumberOfItemsGet(interfaceVarList,numberOfVariables,err,error,*999)
                    DO variableIdx=1,numberOfVariables
                      CALL List_ItemGet(interfaceVarList,variableIdx,variableItem,err,error,*999)
                      interfaceConditionIdx2=variableItem(1)
                      varInterfaceCondition=>solverMapping%interfaceConditions(interfaceConditionIdx2)%ptr
                      IF(ASSOCIATED(varInterfaceCondition)) THEN
                        SELECT CASE(varInterfaceCondition%method)
                        CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD,INTERFACE_CONDITION_PENALTY_METHOD)
                          IF(ASSOCIATED(interfaceCondition%lagrange)) THEN
                            varLagrangeField=>varInterfaceCondition%lagrange%LAGRANGE_FIELD
                            IF(ASSOCIATED(varLagrangeField)) THEN
                              IF(ASSOCIATED(lagrangeField,varLagrangeField)) THEN
                                IF(variableType==variableItem(2)) THEN
                                  variableFound=.TRUE.
                                  EXIT
                                ENDIF
                              ENDIF
                            ELSE
                              CALL FlagError("Variable Lagrange field is not associated.",err,error,*999)
                            ENDIF
                          ELSE
                            CALL FlagError("Variable interface Lagrange is not associated.",err,error,*999)
                          ENDIF
                        CASE(INTERFACE_CONDITION_AUGMENTED_LAGRANGE_METHOD)
                          CALL FlagError("Not implemented.",err,error,*999)
                        CASE(INTERFACE_CONDITION_POINT_TO_POINT_METHOD)
                          CALL FlagError("Not implemented.",err,error,*999)
                        CASE DEFAULT
                          localError="The interface condition method of "// &
                            & TRIM(NumberToVString(varInterfaceCondition%method,"*",err,error))// &
                            & " is invalid."
                          CALL FlagError(localError,err,error,*999)
                        END SELECT
                      ELSE
                        CALL FlagError("Variable equations set is not associated.",err,error,*999)
                      ENDIF
                    ENDDO !variableIdx
                    IF(.NOT.variableFound) THEN
                      variableItem(1)=interfaceConditionIdx
                      variableItem(2)=variableType
                      CALL List_ItemAdd(interfaceVarList,variableItem,err,error,*999)
                    ENDIF
                  ENDIF
                ELSE
                  CALL FlagError("Lagrange field is not associated.",err,error,*999)
                ENDIF
              ELSE
                CALL FlagError("Interface condition Lagrange is not asssociated.",err,error,*999)
              ENDIF
            CASE(INTERFACE_CONDITION_AUGMENTED_LAGRANGE_METHOD)
              CALL FlagError("Not implemented.",err,error,*999)
            CASE(INTERFACE_CONDITION_POINT_TO_POINT_METHOD)
              CALL FlagError("Not implemented.",err,error,*999)
            CASE DEFAULT
              localError="The interface condition method of "// &
                & TRIM(NumberToVString(interfaceCondition%method,"*",err,error))// &
                & " is invalid."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          ELSE
            CALL FlagError("Interface condition is not associated.",err,error,*999)
          ENDIF
        ELSE
          localError="The specified interface condition index of "// &
            & TRIM(NumberToVString(interfaceConditionIdx,"*",err,error))// &
            & " is invalid. The index must be > 0 and <= "// &
            & TRIM(NumberToVString(solverMapping%numberOfInterfaceConditions,"*",err,error))//"."        
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Interface variable list is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Solver mapping is not associated.",err,error,*999)
    ENDIF
    
    CALL Exits("SolverMapping_CreateValuesCacheInterfaceVarListAdd")
    RETURN
999 CALL Errors("SolverMapping_CreateValuesCacheInterfaceVarListAdd",err,error)    
    CALL Exits("SolverMapping_CreateValuesCacheInterfaceVarListAdd")
    RETURN 1
   
  END SUBROUTINE SolverMapping_CreateValuesCacheInterfaceVarListAdd

  !
  !================================================================================================================================
  !

  !>Destroy a solver mapping.
  SUBROUTINE SolverMapping_Destroy(solverMapping,err,error,*)

    !Argument variables
    TYPE(SolverMappingType), POINTER :: solverMapping !<A pointer the solver mapping to destroy
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    CALL Enters("SolverMapping_Destroy",err,error,*999)

    IF(ASSOCIATED(solverMapping)) THEN
      CALL SolverMapping_Finalise(solverMapping,err,error,*999)
    ELSE
      CALL FlagError("Solver mapping is not associated.",err,error,*999)
    ENDIF
        
    CALL Exits("SolverMapping_Destroy")
    RETURN
999 CALL Errors("SolverMapping_Destroy",err,error)    
    CALL Exits("SolverMapping_Destroy")
    RETURN 1
   
  END SUBROUTINE SolverMapping_Destroy

  !
  !================================================================================================================================
  !

  !>Finalises a equations column to solver columns map and deallocates all memory.
  SUBROUTINE SolverMapping_EquationsColToSolverColsMapFinalise(equationsColToSolverColsMap,err,error,*)

    !Argument variables
    TYPE(EquationsColToSolverColsMapType), INTENT(INOUT) :: equationsColToSolverColsMap !<The equations col to solver cols map to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    CALL Enters("SolverMapping_EquationsColToSolverColsMapFinalise",err,error,*999)

    IF(ALLOCATED(equationsColToSolverColsMap%solverCols)) DEALLOCATE(equationsColToSolverColsMap%solverCols)
    IF(ALLOCATED(equationsColToSolverColsMap%couplingCoefficients)) DEALLOCATE(equationsColToSolverColsMap%couplingCoefficients)
        
    CALL Exits("SolverMapping_EquationsColToSolverColsMapFinalise")
    RETURN
999 CALL Errors("SolverMapping_EquationsColToSolverColsMapFinalise",err,error)    
    CALL Exits("SolverMapping_EquationsColToSolverColsMapFinalise")
    RETURN 1
   
  END SUBROUTINE SolverMapping_EquationsColToSolverColsMapFinalise

  !
  !================================================================================================================================
  !

  !>Initialises an equations column to solver columns map
  SUBROUTINE SolverMapping_EquationsColToSolverColsMapInitialise(equationsColToSolverColsMap,err,error,*)

    !Argument variables
    TYPE(EquationsColToSolverColsMapType), INTENT(OUT) :: equationsColToSolverColsMap !<The equations column to solver columns map to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    CALL Enters("SolverMapping_EquationsColToSolverColsMapInitialise",err,error,*999)

    equationsColToSolverColsMap%numberOfSolverCols=0
    
    CALL Exits("SolverMapping_EquationsColToSolverColsMapInitialise")
    RETURN
999 CALL Errors("SolverMapping_EquationsColToSolverColsMapInitialise",err,error)    
    CALL Exits("SolverMapping_EquationsColToSolverColsMapInitialise")
    RETURN 1
   
  END SUBROUTINE SolverMapping_EquationsColToSolverColsMapInitialise

  !
  !================================================================================================================================
  !

  !>Sets/changes the mapping of global variables to a solver matrix for the solver mapping
  SUBROUTINE SolverMapping_EquationsVariablesToSolverMatrixSet0(solverMapping,solverMatrixIndex,equationsSetIndex,variableType, &
    & err,error,*)

    !Argument variables
    TYPE(SolverMappingType), POINTER :: solverMapping !<A pointer to the solver mapping
    INTEGER(INTG), INTENT(IN) :: solverMatrixIndex !<The solver matrix Index to set the equations variables for
    INTEGER(INTG), INTENT(IN) :: equationsSetIndex !<The equations set index in the solver mapping to specify the variable types for
    INTEGER(INTG), INTENT(IN) :: variableType !<The variable types to map to the solver matrix
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    
    CALL Enters("SolverMapping_EquationsVariablesToSolverMatrixSet0",err,error,*999)

    CALL SolverMapping_EquationsVariablesToSolverMatrixSet1(solverMapping,solverMatrixIndex,equationsSetIndex,[variableType], &
      & err,error,*999)

    CALL Exits("SolverMapping_EquationsVariablesToSolverMatrixSet0")
    RETURN
999 CALL Errors("SolverMapping_EquationsVariablesToSolverMatrixSet0",err,error)
    CALL Exits("SolverMapping_EquationsVariablesToSolverMatrixSet0")
    RETURN 1
  END SUBROUTINE SolverMapping_EquationsVariablesToSolverMatrixSet0
    
  !
  !================================================================================================================================
  !

  !>Sets/changes the mapping of global variables to a solver matrix for the solver mapping
  SUBROUTINE SolverMapping_EquationsVariablesToSolverMatrixSet1(solverMapping,solverMatrixIndex,equationsSetIndex,variableTypes, &
    & err,error,*)

    !Argument variables
    TYPE(SolverMappingType), POINTER :: solverMapping !<A pointer to the solver mapping
    INTEGER(INTG), INTENT(IN) :: solverMatrixIndex !<The solver matrix Index to set the equations variables for
    INTEGER(INTG), INTENT(IN) :: equationsSetIndex !<The equations set index in the solver mapping to specify the variable types for
    INTEGER(INTG), INTENT(IN) :: variableTypes(:) !<The variable types to map to the solver matrix
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: variableIdx,variableType
    TYPE(EQUATIONS_TYPE), POINTER :: equations
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: equationsMapping
    TYPE(EQUATIONS_MAPPING_LINEAR_TYPE), POINTER :: linearMapping
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet
    TYPE(FIELD_TYPE), POINTER :: dependentField
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: dependentVariable
    TYPE(VARYING_STRING) :: localError
    
    CALL Enters("SolverMapping_EquationsVariablesToSolverMatrixSet1",err,error,*999)

    IF(ASSOCIATED(solverMapping)) THEN
      IF(solverMapping%solverMappingFinished) THEN
        CALL FlagError("Solver mappings has been finished.",err,error,*999)
      ELSE
        IF(ASSOCIATED(solverMapping%createValuesCache)) THEN
          IF(solverMatrixIndex>=1.AND.solverMatrixIndex<=solverMapping%numberOfSolverMatrices) THEN
            IF(equationsSetIndex>=1.AND.equationsSetIndex<=solverMapping%numberOfEquationsSets) THEN
              equationsSet=>solverMapping%equationsSets(equationsSetIndex)%ptr
              IF(ASSOCIATED(equationsSet)) THEN
                dependentField=>equationsSet%dependent%DEPENDENT_FIELD
                IF(ASSOCIATED(dependentField)) THEN
                  equations=>equationsSet%equations
                  IF(ASSOCIATED(equations)) THEN
                    equationsMapping=>equations%EQUATIONS_MAPPING
                    IF(ASSOCIATED(equationsMapping)) THEN
                      linearMapping=>equationsMapping%LINEAR_MAPPING
                      IF(ASSOCIATED(linearMapping)) THEN
                        IF(SIZE(variableTypes,1)>=1.AND.SIZE(variableTypes,1)<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
                          DO variableIdx=1,SIZE(variableTypes,1)
!!\todo CHECK THAT THE VARIABLE TYPE IS NOT REPEATED
                            variableType=variableTypes(variableIdx)
                            CALL Field_VariableGet(dependentField,variableType,dependentVariable,err,error,*999)
                            IF(linearMapping%VAR_TO_EQUATIONS_MATRICES_MAPS(variableType)%NUMBER_OF_EQUATIONS_MATRICES==0) THEN
                              localError="The variable type of "// &
                                & TRIM(NumberToVString(variableType,"*",err,error))// &
                                & " at position "//TRIM(NumberToVString(variableIdx,"*",err,error))// &
                                & " in the array is invalid. That variable type is not mapped to any equations matrices."
                            ENDIF
                          ENDDO !variableIdx
                          solverMapping%createValuesCache%matrixVariableTypes(0,equationsSetIndex,solverMatrixIndex)= &
                            & SIZE(variableTypes,1)
                          solverMapping%createValuesCache%matrixVariableTypes(1:SIZE(variableTypes,1),equationsSetIndex, &
                            & solverMatrixIndex)=variableTypes
                        ELSE
                          localError="The supplied size of variable types array of "// &
                            & TRIM(NumberToVString(SIZE(variableTypes,1),"*",err,error))// &
                            & " is invalid. The size must be between 1 and "// &
                            & TRIM(NumberToVString(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",err,error))//"."
                          CALL FlagError(localError,err,error,*999)
                        ENDIF
                      ELSE
                        CALL FlagError("Equations mapping linear mapping is not associated.",err,error,*999)
                      ENDIF
                    ELSE
                      CALL FlagError("Equations mapping is not associated.",err,error,*999)
                    ENDIF
                  ELSE
                    CALL FlagError("Equations is not associated.",err,error,*999)
                  ENDIF
                ELSE
                  CALL FlagError("Equations set dependent field is not associated.",err,error,*999)
                ENDIF
              ELSE
                CALL FlagError("Equations set is not associated.",err,error,*999)
              ENDIF
            ELSE
              localError="The equations set index of "//TRIM(NumberToVString(equationsSetIndex,"*",err,error))// &
                & " is invalid. The number must be >= 1 and <= "// &
                & TRIM(NumberToVString(solverMapping%numberOfEquationsSets,"*",err,error))//"."
              CALL FlagError(localError,err,error,*999)
            ENDIF
          ELSE
            localError="The solver matrix number of "//TRIM(NumberToVString(solverMatrixIndex,"*",err,error))// &
              & " is invalid. The number must be >= 1 and <= "// &
              & TRIM(NumberToVString(solverMapping%numberOfSolverMatrices,"*",err,error))//"."
            CALL FlagError(localError,err,error,*999)
          ENDIF
        ELSE
          CALL FlagError("Solver mapping create values cache is not associated",err,error,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Solver mapping is not associated",err,error,*999)
    ENDIF
    
    CALL Exits("SolverMapping_EquationsVariablesToSolverMatrixSet1")
    RETURN
999 CALL Errors("SolverMapping_EquationsVariablesToSolverMatrixSet1",err,error)
    CALL Exits("SolverMapping_EquationsVariablesToSolverMatrixSet1")
    RETURN 1
  END SUBROUTINE SolverMapping_EquationsVariablesToSolverMatrixSet1
  
  !
  !================================================================================================================================
  !

  !>Finalises a equations row to solver rows map and deallocates all memory.
  SUBROUTINE SolverMapping_EquationsRowToSolverRowsMapFinalise(equationsRowToSolverRowsMap,err,error,*)

    !Argument variables
    TYPE(EquationsRowToSolverRowsMapType), INTENT(INOUT) :: equationsRowToSolverRowsMap !<The equations row to solver rows map to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    CALL Enters("SolverMapping_EquationsRowToSolverRowsMapFinalise",err,error,*999)

    IF(ALLOCATED(equationsRowToSolverRowsMap%solverRows)) DEALLOCATE(equationsRowToSolverRowsMap%solverRows)
    IF(ALLOCATED(equationsRowToSolverRowsMap%couplingCoefficients)) DEALLOCATE(equationsRowToSolverRowsMap%couplingCoefficients)
        
    CALL Exits("SolverMapping_EquationsRowToSolverRowsMapFinalise")
    RETURN
999 CALL Errors("SolverMapping_EquationsRowToSolverRowsMapFinalise",err,error)    
    CALL Exits("SolverMapping_EquationsRowToSolverRowsMapFinalise")
    RETURN 1
   
  END SUBROUTINE SolverMapping_EquationsRowToSolverRowsMapFinalise

  !
  !================================================================================================================================
  !

  !>Initialises an equations row to solver rows map
  SUBROUTINE SolverMapping_EquationsRowToSolverRowsMapInitialise(equationsRowToSolverRowsMap,err,error,*)

    !Argument variables
    TYPE(EquationsRowToSolverRowsMapType), INTENT(OUT) :: equationsRowToSolverRowsMap !<The equations row to solver rows map to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    CALL Enters("SolverMapping_EquationsRowToSolverRowsMapInitialise",err,error,*999)

    equationsRowToSolverRowsMap%numberOfSolverRows=0
    
    CALL Exits("SolverMapping_EquationsRowToSolverRowsMapInitialise")
    RETURN
999 CALL Errors("SolverMapping_EquationsRowToSolverRowsMapInitialise",err,error)    
    CALL Exits("SolverMapping_EquationsRowToSolverRowsMapInitialise")
    RETURN 1
   
  END SUBROUTINE SolverMapping_EquationsRowToSolverRowsMapInitialise

  !
  !================================================================================================================================
  !

  !>Adds an equations set to a solver mapping
  SUBROUTINE SolverMapping_EquationsSetAdd(solverMapping,equationsSet,equationsSetIndex,err,error,*)

    !Argument variables
    TYPE(SolverMappingType), POINTER :: solverMapping !<A pointer the solver mapping to add the equations set to
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set to add
    INTEGER(INTG), INTENT(OUT) :: equationsSetIndex !<On exit, the index of the equations set in the solver mapping
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: equationsSetIdx,matrixIdx,solverMatrixIdx,variableIdx,variableType
    INTEGER(INTG), ALLOCATABLE :: newDynamicVariableType(:),newMatrixVariableTypes(:,:,:),newLhsVariableType(:), &
      & newRhsVariableType(:),newResidualVariableTypes(:,:),newSourceVariableType(:)
    LOGICAL :: matrixDone
    TYPE(EQUATIONS_TYPE), POINTER :: equations
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: equationsMapping
    TYPE(EQUATIONS_SET_ptr_TYPE), ALLOCATABLE :: newEquationsSets(:)
    TYPE(FIELD_TYPE), POINTER :: dependentField
    TYPE(LIST_ptr_TYPE), POINTER :: newInterfaceIndices(:)
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: solverEquations
    TYPE(VARYING_STRING) :: localError

    NULLIFY(newInterfaceIndices)
   
    CALL Enters("SolverMapping_EquationsSetAdd",err,error,*999)

    equationsSetIdx=0
    IF(ASSOCIATED(solverMapping)) THEN
      IF(solverMapping%solverMappingFinished) THEN
        CALL FlagError("Solver mapping has been finished.",err,error,*999)
      ELSE
        solverEquations=>solverMapping%solverEquations
        IF(ASSOCIATED(solverEquations)) THEN
          IF(ASSOCIATED(equationsSet)) THEN
            IF(equationsSet%EQUATIONS_SET_FINISHED) THEN
              equations=>equationsSet%equations
              IF(ASSOCIATED(equations)) THEN
                equationsMapping=>equations%EQUATIONS_MAPPING
                IF(ASSOCIATED(equationsMapping)) THEN                
                  IF(ASSOCIATED(solverMapping%createValuesCache)) THEN
                    IF(solverMapping%numberOfEquationsSets>0) THEN
                      ALLOCATE(newDynamicVariableType(solverMapping%numberOfEquationsSets+1),STAT=err)
                      IF(err/=0) CALL FlagError("Could not allocate new dynamic variable type.",err,error,*999)
                      ALLOCATE(newMatrixVariableTypes(0:FIELD_NUMBER_OF_VARIABLE_TYPES,solverMapping% &
                        & numberOfEquationsSets+1,solverMapping%numberOfSolverMatrices),STAT=err)
                      IF(err/=0) CALL FlagError("Could not allocate new matrix variable types.",err,error,*999)
                      ALLOCATE(newResidualVariableTypes(0:FIELD_NUMBER_OF_VARIABLE_TYPES, &
                        & solverMapping%numberOfEquationsSets+1),STAT=err)
                      IF(err/=0) CALL FlagError("Could not allocate new residual variable type.",err,error,*999)
                      ALLOCATE(newLhsVariableType(solverMapping%numberOfEquationsSets+1),STAT=err)
                      IF(err/=0) CALL FlagError("Could not allocate new LHS variable type.",err,error,*999)
                      ALLOCATE(newRhsVariableType(solverMapping%numberOfEquationsSets+1),STAT=err)
                      IF(err/=0) CALL FlagError("Could not allocate new RHS variable type.",err,error,*999)
                      ALLOCATE(newSourceVariableType(solverMapping%numberOfEquationsSets+1),STAT=err)
                      IF(err/=0) CALL FlagError("Could not allocate new source variable type.",err,error,*999)
                      ALLOCATE(newInterfaceIndices(solverMapping%numberOfEquationsSets+1),STAT=err)
                      IF(err/=0) CALL FlagError("Could not allocate new interface indices.",err,error,*999)
                      DO equationsSetIdx=1,solverMapping%numberOfEquationsSets
                        newInterfaceIndices(equationsSetIdx)%ptr=>solverMapping%createValuesCache%interfaceIndices( &
                          & equationsSetIdx)%ptr
                      ENDDO !equationsSetIdx
                      ALLOCATE(newEquationsSets(solverMapping%numberOfEquationsSets+1),STAT=err)
                      IF(err/=0) CALL FlagError("Could not allocate new equations sets.",err,error,*999)

                      IF(ASSOCIATED(solverMapping%createValuesCache%interfaceIndices)) &
                        & DEALLOCATE(solverMapping%createValuesCache%interfaceIndices)

                      newDynamicVariableType(1:solverMapping%numberOfEquationsSets)= &
                        & solverMapping%createValuesCache%dynamicVariableType
                      newMatrixVariableTypes(:,1:solverMapping%numberOfEquationsSets,:)= &
                        & solverMapping%createValuesCache%matrixVariableTypes
                      newResidualVariableTypes(:,1:solverMapping%numberOfEquationsSets)= &
                        & solverMapping%createValuesCache%residualVariableTypes
                      newLhsVariableType(1:solverMapping%numberOfEquationsSets)= &
                        & solverMapping%createValuesCache%lhsVariableType
                      newRhsVariableType(1:solverMapping%numberOfEquationsSets)= &
                        & solverMapping%createValuesCache%rhsVariableType
                      newSourceVariableType(1:solverMapping%numberOfEquationsSets)= &
                        & solverMapping%createValuesCache%sourceVariableType
                      DO equationsSetIdx=1,solverMapping%numberOfEquationsSets
                        newEquationsSets(equationsSetIdx)%ptr=>solverMapping%equationsSets(equationsSetIdx)%ptr
                      ENDDO !equationsSetIdx
                      solverMapping%createValuesCache%interfaceIndices=>newInterfaceIndices

                      newDynamicVariableType(solverMapping%numberOfEquationsSets+1)=0
                      newMatrixVariableTypes(:,solverMapping%numberOfEquationsSets+1,:)=0
                      newResidualVariableTypes(:,solverMapping%numberOfEquationsSets+1)=0
                      newLhsVariableType(solverMapping%numberOfEquationsSets+1)=0
                      newRhsVariableType(solverMapping%numberOfEquationsSets+1)=0
                      newSourceVariableType(solverMapping%numberOfEquationsSets+1)=0

                      CALL MOVE_ALLOC(newDynamicVariableType,solverMapping%createValuesCache%dynamicVariableType)
                      CALL MOVE_ALLOC(newMatrixVariableTypes,solverMapping%createValuesCache%matrixVariableTypes)
                      CALL MOVE_ALLOC(newResidualVariableTypes,solverMapping%createValuesCache%residualVariableTypes)
                      CALL MOVE_ALLOC(newLhsVariableType,solverMapping%createValuesCache%lhsVariableType)
                      CALL MOVE_ALLOC(newRhsVariableType,solverMapping%createValuesCache%rhsVariableType)
                      CALL MOVE_ALLOC(newSourceVariableType,solverMapping%createValuesCache%sourceVariableType)
                      CALL MOVE_ALLOC(newEquationsSets,solverMapping%equationsSets)
                      solverMapping%createValuesCache%interfaceIndices=>newInterfaceIndices

                    ELSE IF(solverMapping%numberOfEquationsSets==0) THEN
                      
                      ALLOCATE(newDynamicVariableType(solverMapping%numberOfEquationsSets+1),STAT=err)
                      IF(err/=0) CALL FlagError("Could not allocate dynamic variable type.",err,error,*999)
                      ALLOCATE(newMatrixVariableTypes(0:FIELD_NUMBER_OF_VARIABLE_TYPES, &
                        & solverMapping%numberOfEquationsSets+1,solverMapping%numberOfSolverMatrices),STAT=err)
                      IF(err/=0) CALL FlagError("Could not allocate matrix variable types.",err,error,*999)
                      ALLOCATE(newResidualVariableTypes(0:FIELD_NUMBER_OF_VARIABLE_TYPES, &
                        & solverMapping%numberOfEquationsSets+1),STAT=err)
                      IF(err/=0) CALL FlagError("Could not allocate residual variable type.",err,error,*999)
                      ALLOCATE(newLhsVariableType(solverMapping%numberOfEquationsSets+1),STAT=err)
                      IF(err/=0) CALL FlagError("Could not allocate LHS variable type.",err,error,*999)
                      ALLOCATE(newRhsVariableType(solverMapping%numberOfEquationsSets+1),STAT=err)
                      IF(err/=0) CALL FlagError("Could not allocate RHS variable type.",err,error,*999)
                      ALLOCATE(newSourceVariableType(solverMapping%numberOfEquationsSets+1),STAT=err)
                      IF(err/=0) CALL FlagError("Could not allocate source variable type.",err,error,*999)
                      ALLOCATE(newInterfaceIndices(solverMapping%numberOfEquationsSets+1),STAT=err)
                      IF(err/=0) CALL FlagError("Could not allocate new interface indices.",err,error,*999)
                      ALLOCATE(newEquationsSets(solverMapping%numberOfEquationsSets+1),STAT=err)
                      IF(err/=0) CALL FlagError("Could not allocate new equations sets.",err,error,*999)

                      newDynamicVariableType=0
                      newMatrixVariableTypes=0
                      newResidualVariableTypes=0
                      newLhsVariableType=0
                      newRhsVariableType=0
                      newSourceVariableType=0

                      CALL MOVE_ALLOC(newDynamicVariableType,solverMapping%createValuesCache%dynamicVariableType)
                      CALL MOVE_ALLOC(newMatrixVariableTypes,solverMapping%createValuesCache%matrixVariableTypes)
                      CALL MOVE_ALLOC(newResidualVariableTypes,solverMapping%createValuesCache%residualVariableTypes)
                      CALL MOVE_ALLOC(newLhsVariableType,solverMapping%createValuesCache%lhsVariableType)
                      CALL MOVE_ALLOC(newRhsVariableType,solverMapping%createValuesCache%rhsVariableType)
                      CALL MOVE_ALLOC(newSourceVariableType,solverMapping%createValuesCache%sourceVariableType)
                      CALL MOVE_ALLOC(newEquationsSets,solverMapping%equationsSets)
                      IF(ASSOCIATED(solverMapping%createValuesCache%interfaceIndices)) &
                        & DEALLOCATE(solverMapping%createValuesCache%interfaceIndices)
                      solverMapping%createValuesCache%interfaceIndices=>newInterfaceIndices
                    ELSE
                      CALL FlagError("The number of equations sets is < 0.",err,error,*999)
                    ENDIF
                    NULLIFY(solverMapping%createValuesCache%interfaceIndices(solverMapping%numberOfEquationsSets+1)%ptr)
                    CALL List_CreateStart(solverMapping%createValuesCache%interfaceIndices(solverMapping% &
                      & numberOfEquationsSets+1)%ptr,err,error,*999)
                    CALL List_DataTypeSet(solverMapping%createValuesCache%interfaceIndices(solverMapping% &
                       & numberOfEquationsSets+1)%ptr,LIST_INTG_TYPE,err,error,*999)
                    CALL List_DataDimensionSet(solverMapping%createValuesCache%interfaceIndices(solverMapping% &
                       & numberOfEquationsSets+1)%ptr,2,err,error,*999)
                    CALL List_KeyDimensionSet(solverMapping%createValuesCache%interfaceIndices(solverMapping% &
                      & numberOfEquationsSets+1)%ptr,1,err,error,*999)
                    CALL List_CreateFinish(solverMapping%createValuesCache%interfaceIndices(solverMapping% &
                      & numberOfEquationsSets+1)%ptr,err,error,*999)
                    SELECT CASE(equations%TIME_DEPENDENCE)
                    CASE(EQUATIONS_STATIC,EQUATIONS_QUASISTATIC)
                      SELECT CASE(equations%linearity)
                      CASE(EQUATIONS_LINEAR)
                        IF(ASSOCIATED(equationsMapping%LINEAR_MAPPING)) THEN
                          !Linear matrices to map. 
                          !Map the first matrix variable found in the equations set to the first solver matrix, the second
                          !variable found to the second, etc.
                          variableType=1
                          DO matrixIdx=1,solverMapping%numberOfSolverMatrices
                            matrixDone=.FALSE.
                            DO WHILE(variableType<=FIELD_NUMBER_OF_VARIABLE_TYPES.AND..NOT.matrixDone)
                              IF(equationsMapping%LINEAR_MAPPING%VAR_TO_EQUATIONS_MATRICES_MAPS(variableType)% &
                                & NUMBER_OF_EQUATIONS_MATRICES>0) THEN                  
                                solverMapping%createValuesCache%matrixVariableTypes(0, &
                                  & solverMapping%numberOfEquationsSets+1,matrixIdx)=1
                                solverMapping%createValuesCache%matrixVariableTypes(1, &
                                  & solverMapping%numberOfEquationsSets+1,matrixIdx)=variableType
                                matrixDone=.TRUE.
                              ELSE
                                variableType=variableType+1
                              ENDIF
                            ENDDO
                            IF(.NOT.matrixDone) THEN
                              !Error - could not find any more variables to map to this solver matrix
                              localError="Could not find any unmapped variables for solver matrix "// &
                                & TRIM(NumberToVString(matrixIdx,"*",err,error))//"."
                              CALL FlagError(localError,err,error,*999)
                            ENDIF
                          ENDDO !matrixIdx
                          !Check if there are still unmapped matrix variables.
                          DO variableIdx=variableType+1,FIELD_NUMBER_OF_VARIABLE_TYPES
                            IF(equationsMapping%LINEAR_MAPPING%VAR_TO_EQUATIONS_MATRICES_MAPS(variableIdx)% &
                              & NUMBER_OF_EQUATIONS_MATRICES>0) THEN
                              localError="Variable type "//TRIM(NumberToVString(variableIdx,"*",err,error))// &
                                & " is mapped to a linear matrix but has not been mapped to any solver matrices."
                              CALL FlagError(localError,err,error,*999)
                            ENDIF
                          ENDDO !variableIdx
                        ELSE
                          CALL FlagError("Equations mapping linear mapping is not associated.",err,error,*999)
                        ENDIF
                      CASE(EQUATIONS_NONLINEAR)
                        IF(ASSOCIATED(equationsMapping%NONLINEAR_MAPPING)) THEN
                          !Set the number of residual variables for this equations set
                          solverMapping%createValuesCache%residualVariableTypes(0,solverMapping%numberOfEquationsSets+1)= &
                            & equationsMapping%NONLINEAR_MAPPING%NUMBER_OF_RESIDUAL_VARIABLES
                          !Map the residual variables to the solver Jacobian
                          DO matrixIdx=1,equationsMapping%NONLINEAR_MAPPING%NUMBER_OF_RESIDUAL_VARIABLES
                            solverMapping%createValuesCache%residualVariableTypes(matrixIdx,solverMapping% &
                              & numberOfEquationsSets+1)=equationsMapping%NONLINEAR_MAPPING% &
                              & JACOBIAN_TO_VAR_MAP(matrixIdx)%VARIABLE_TYPE
                          ENDDO
                          IF(ASSOCIATED(equationsMapping%LINEAR_MAPPING)) THEN
                            !If there are linear matrices operating on the residual variable then map them to the
                            !solver matrix (Jacobian)
                            IF(solverMapping%numberOfSolverMatrices==1) THEN
                              DO matrixIdx=1,equationsMapping%NONLINEAR_MAPPING%NUMBER_OF_RESIDUAL_VARIABLES
                                variableType=equationsMapping%NONLINEAR_MAPPING%JACOBIAN_TO_VAR_MAP(matrixIdx)%VARIABLE_TYPE
                                IF(equationsMapping%LINEAR_MAPPING%VAR_TO_EQUATIONS_MATRICES_MAPS(variableType)% &
                                  & NUMBER_OF_EQUATIONS_MATRICES>0) THEN
                                  solverMapping%createValuesCache%matrixVariableTypes(0, &
                                    & solvermapping%numberOfEquationsSets+1,1)=1
                                  solverMapping%createValuesCache%matrixVariableTypes(1, &
                                    & solverMapping%numberOfEquationsSets+1,1)=variableType
                                ENDIF
                              ENDDO !matrixIdx
                            ELSE
                              localError="Invalid number of solve matrices. For nonlinear solver equations there should "// &
                                & "be 1 solver matrix and there are "// &
                                & TRIM(NumberToVString(solverMapping%numberOfSolverMatrices,"*",err,error))// &
                                & " solver matrices."
                              CALL FlagError(localError,err,error,*999)
                            ENDIF
                          ENDIF
                        ELSE
                          CALL FlagError("Equations mapping nonlinear mapping is not associated.",err,error,*999)
                        ENDIF
                      CASE DEFAULT
                        localError="The equations linearity type of "// &
                          & TRIM(NumberToVString(solverEquations%linearity,"*",err,error))//" is invalid."
                        CALL FlagError(localError,err,error,*999)
                      END SELECT
                    CASE(EQUATIONS_FIRST_ORDER_DYNAMIC,EQUATIONS_SECOND_ORDER_DYNAMIC)
                      SELECT CASE(equations%linearity)
                      CASE(EQUATIONS_LINEAR)
                        IF(ASSOCIATED(equationsMapping%DYNAMIC_MAPPING)) THEN
                          solverMapping%createValuesCache%dynamicVariableType(solverMapping%numberOfEquationsSets+1)= &
                            & equationsMapping%DYNAMIC_MAPPING%DYNAMIC_VARIABLE_TYPE
                        ELSE
                          CALL FlagError("Equations mapping dynamic mapping is not associated.",err,error,*999)
                        ENDIF
                      CASE(EQUATIONS_NONLINEAR)
                        IF(ASSOCIATED(equationsMapping%DYNAMIC_MAPPING)) THEN
                          solverMapping%createValuesCache%dynamicVariableType(solverMapping%numberOfEquationsSets+1)= &
                            & equationsMapping%DYNAMIC_MAPPING%DYNAMIC_VARIABLE_TYPE
                          DO matrixIdx=1,equationsMapping%NONLINEAR_MAPPING%NUMBER_OF_RESIDUAL_VARIABLES
                            solverMapping%createValuesCache%residualVariableTypes(matrixIdx,solverMapping% &
                              & numberOfEquationsSets+1)=equationsMapping%NONLINEAR_MAPPING% &
                              & JACOBIAN_TO_VAR_MAP(matrixIdx)%VARIABLE_TYPE
                          ENDDO !matrixIdx
                        ELSE
                          CALL FlagError("Equations mapping dynamic mapping is not associated.",err,error,*999)
                        ENDIF
                      CASE DEFAULT
                        localError="The equations linearity type of "// &
                          & TRIM(NumberToVString(solverEquations%linearity,"*",err,error))//" is invalid."
                        CALL FlagError(localError,err,error,*999)
                      END SELECT
                    CASE DEFAULT
                      localError="The equations time dependence type of "// &
                        & TRIM(NumberToVString(solverEquations%TIME_DEPENDENCE,"*",err,error))//" is invalid."
                      CALL FlagError(localError,err,error,*999)
                    END SELECT
                    IF(ASSOCIATED(equationsMapping%lhsMapping)) THEN
                      solverMapping%createValuesCache%lhsVariableType(solverMapping%numberOfEquationsSets+1)= &
                        & equationsMapping%lhsMapping%lhsVariableType
                    ELSE
                      CALL FlagError("Equations mapping has no LHS mapping.",err,error,*999)
                    ENDIF
                    IF(ASSOCIATED(equationsMapping%RHS_MAPPING)) THEN
                      solverMapping%createValuesCache%rhsVariableType(solverMapping%numberOfEquationsSets+1)= &
                        & equationsMapping%RHS_MAPPING%RHS_VARIABLE_TYPE
                    ENDIF
                    IF(ASSOCIATED(equationsMapping%SOURCE_MAPPING)) THEN
                      solverMapping%createValuesCache%sourceVariableType(solverMapping%numberOfEquationsSets+1)= &
                        & equationsMapping%SOURCE_MAPPING%SOURCE_VARIABLE_TYPE
                    ENDIF
                    solverMapping%equationsSets(solverMapping%numberOfEquationsSets+1)%ptr=>equationsSet
                    solverMapping%numberOfEquationsSets=solverMapping%numberOfEquationsSets+1
                    equationsSetIndex=solverMapping%numberOfEquationsSets
                    
                    !Add the variables to the list of variables
                    dependentField=>equationsSet%dependent%DEPENDENT_FIELD
                    variableType=solverMapping%createValuesCache%dynamicVariableType(equationsSetIndex)
                    CALL SolverMapping_CreateValuesCacheEquationVarListAdd(solverMapping,solverMapping%createValuesCache% &
                      & equationsColVariablesList(1)%ptr,equationsSetIndex,variableType,err,error,*999)
                    DO solverMatrixIdx=1,solverMapping%numberOfSolverMatrices
                      DO matrixIdx=1,solverMapping%createValuesCache%matrixVariableTypes(0,equationsSetIndex,solverMatrixIdx)
                        variableType=solverMapping%createValuesCache%matrixVariableTypes(matrixIdx,equationsSetIndex, &
                          & solverMatrixIdx)
                        CALL SolverMapping_CreateValuesCacheEquationVarListAdd(solverMapping,solverMapping%createValuesCache% &
                          & equationsColVariablesList(1)%ptr,equationsSetIndex,variableType,err,error,*999)
                      ENDDO !matrixIdx
                    ENDDO !solverMatrixIdx
                    DO matrixIdx=1,solverMapping%createValuesCache%residualVariableTypes(0,equationsSetIndex)
                      variableType=solverMapping%createValuesCache%residualVariableTypes(matrixIdx,equationsSetIndex)
                      CALL SolverMapping_CreateValuesCacheEquationVarListAdd(solverMapping,solverMapping%createValuesCache% &
                        & equationsColVariablesList(1)%ptr,equationsSetIndex,variableType,err,error,*999)
                    ENDDO !matrixIdx
                  ELSE
                    CALL FlagError("Solvers mapping create values cache is not associated.",err,error,*999)
                  ENDIF
                ELSE
                  CALL FlagError("Equations mapping is not associated.",err,error,*999)
                ENDIF
              ELSE
                CALL FlagError("Equations set equations is not associated.",err,error,*999)
              ENDIF
            ELSE
              CALL FlagError("Equations set has not been finished.",err,error,*999)
            ENDIF
          ELSE
            CALL FlagError("Equations set is not associated.",err,error,*999)
          ENDIF
        ELSE
          CALL FlagError("Solver mapping solver is not associated.",err,error,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Solver mapping is not associated.",err,error,*999)
    ENDIF
        
    CALL Exits("SolverMapping_EquationsSetAdd")
    RETURN
999 IF(ALLOCATED(newMatrixVariableTypes)) DEALLOCATE(newMatrixVariableTypes)
    IF(ALLOCATED(newResidualVariableTypes)) DEALLOCATE(newResidualVariableTypes)
    IF(ALLOCATED(newLhsVariableType)) DEALLOCATE(newLhsVariableType)
    IF(ALLOCATED(newRhsVariableType)) DEALLOCATE(newRhsVariableType)
    IF(ALLOCATED(newSourceVariableType)) DEALLOCATE(newSourceVariableType)
    IF(ALLOCATED(newEquationsSets)) DEALLOCATE(newEquationsSets)
    CALL Errors("SolverMapping_EquationsSetAdd",err,error)    
    CALL Exits("SolverMapping_EquationsSetAdd")
    RETURN 1
   
  END SUBROUTINE SolverMapping_EquationsSetAdd

  !
  !================================================================================================================================
  !

  !>Finalises a equations set to solver map and deallocates all memory.
  SUBROUTINE SolverMapping_EquationsSetToSolverMapFinalise(equationsSetToSolverMap,err,error,*)

    !Argument variables
    TYPE(EquationsSetToSolverMapType), INTENT(INOUT) :: equationsSetToSolverMap !<The equations set to solver map to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: equationsMatrixIdx,interfaceConditionIdx,rowIdx,solverMatrixIdx
    
    CALL Enters("SolverMapping_EquationsSetToSolverMapFinalise",err,error,*999)

    IF(ALLOCATED(equationsSetToSolverMap%equationsToSolverMatrixMapsInterface)) THEN
      DO interfaceConditionIdx=1,SIZE(equationsSetToSolverMap%equationsToSolverMatrixMapsInterface,1)
        CALL SolverMapping_EquationsToSolverInterfaceFinalise(equationsSetToSolverMap% &
          & equationsToSolverMatrixMapsInterface(interfaceConditionIdx),err,error,*999)
      ENDDO !interfaceConditionIdx
      DEALLOCATE(equationsSetToSolverMap%equationsToSolverMatrixMapsInterface)
    ENDIF
    IF(ALLOCATED(equationsSetToSolverMap%equationsToSolverMatrixMapsSm)) THEN
      DO solverMatrixIdx=1,SIZE(equationsSetToSolverMap%equationsToSolverMatrixMapsSm,1)
        CALL SolverMapping_EquationsToSolverMatrixMapsSmFinalise(equationsSetToSolverMap% &
          & equationsToSolverMatrixMapsSm(solverMatrixIdx),err,error,*999)
      ENDDO !solverMatrixIdx
      DEALLOCATE(equationsSetToSolverMap%equationsToSolverMatrixMapsSm)
    ENDIF
    IF(ALLOCATED(equationsSetToSolverMap%equationsToSolverMatrixMapsEm)) THEN
      DO equationsMatrixIdx=1,SIZE(equationsSetToSolverMap%equationsToSolverMatrixMapsEm,1)
        CALL SolverMapping_EquationsToSolverMatrixMapsEmFinalise(equationsSetToSolverMap% &
          & equationsToSolverMatrixMapsEm(equationsMatrixIdx),err,error,*999)
      ENDDO !equationsMatrixIdx
      DEALLOCATE(equationsSetToSolverMap%equationsToSolverMatrixMapsEm)
    ENDIF
    IF(ALLOCATED(equationsSetToSolverMap%equationsToSolverMatrixMapsJm)) THEN
      CALL SolverMapping_EquationsToSolverMatrixMapsJmFinalise(equationsSetToSolverMap% &
        & equationsToSolverMatrixMapsJm,err,error,*999)
      DEALLOCATE(equationsSetToSolverMap%equationsToSolverMatrixMapsJm)
    ENDIF
    IF(ALLOCATED(equationsSetToSolverMap%equationsRowToSolverRowsMaps)) THEN
      DO rowIdx=1,SIZE(equationsSetToSolverMap%equationsRowToSolverRowsMaps,1)
        CALL SolverMapping_EquationsRowToSolverRowsMapFinalise(equationsSetToSolverMap% &
          & equationsRowToSolverRowsMaps(rowIdx),err,error,*999)
      ENDDO !rowIdx
      DEALLOCATE(equationsSetToSolverMap%equationsRowToSolverRowsMaps)
    ENDIF
        
    CALL Exits("SolverMapping_EquationsSetToSolverMapFinalise")
    RETURN
999 CALL Errors("SolverMapping_EquationsSetToSolverMapFinalise",err,error)    
    CALL Exits("SolverMapping_EquationsSetToSolverMapFinalise")
    RETURN 1
   
  END SUBROUTINE SolverMapping_EquationsSetToSolverMapFinalise

  !
  !================================================================================================================================
  !

  !>Initialises a equations set to solver map.
  SUBROUTINE SolverMapping_EquationsSetToSolverMapInitialise(equationsSetToSolverMap,err,error,*)

    !Argument variables
    TYPE(EquationsSetToSolverMapType), INTENT(OUT) :: equationsSetToSolverMap !<The equations set to solver map to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    CALL Enters("SolverMapping_EquationsSetToSolverMapInitialise",err,error,*999)

    equationsSetToSolverMap%equationsSetIndex=0
    NULLIFY(equationsSetToSolverMap%solverMapping)
    NULLIFY(equationsSetToSolverMap%equations)
    equationsSetToSolverMap%numberOfInterfaceConditions=0

    CALL Exits("SolverMapping_EquationsSetToSolverMapInitialise")
    RETURN
999 CALL Errors("SolverMapping_EquationsSetToSolverMapInitialise",err,error)    
    CALL Exits("SolverMapping_EquationsSetToSolverMapInitialise")
    RETURN 1
   
  END SUBROUTINE SolverMapping_EquationsSetToSolverMapInitialise

  !
  !================================================================================================================================
  !

  !>Finalises a equations set to solver matrix map and deallocates all memory.
  SUBROUTINE SolverMapping_EquationsToSolverMapsFinalise(equationsToSolverMap,err,error,*)

    !Argument variables
    TYPE(EquationsToSolverMapsType), POINTER :: equationsToSolverMap !<The equations set to solver map to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: columnIdx
    
    CALL Enters("SolverMapping_EquationsToSolverMapsFinalise",err,error,*999)

    IF(ASSOCIATED(equationsToSolverMap)) THEN
      IF(ALLOCATED(equationsToSolverMap%equationsColToSolverColsMap)) THEN
        DO columnIdx=1,SIZE(equationsToSolverMap%equationsColToSolverColsMap,1)
          CALL SolverMapping_EquationsColToSolverColsMapFinalise(equationsToSolverMap% &
            & equationsColToSolverColsMap(columnIdx),err,error,*999)
        ENDDO !columnIdx
        DEALLOCATE(equationsToSolverMap%equationsColToSolverColsMap)
      ENDIF
    ENDIF
        
    CALL Exits("SolverMapping_EquationsToSolverMapsFinalise")
    RETURN
999 CALL Errors("SolverMapping_EquationsToSolverMapsFinalise",err,error)    
    CALL Exits("SolverMapping_EquationsToSolverMapsFinalise")
    RETURN 1
   
  END SUBROUTINE SolverMapping_EquationsToSolverMapsFinalise

  !
  !================================================================================================================================
  !

  !>Initialises an equations to solver maps
  SUBROUTINE SolverMapping_EquationsToSolverMapsInitialise(equationsToSolverMaps,err,error,*)

    !Argument variables
    TYPE(EquationsToSolverMapsType), POINTER :: equationsToSolverMaps !<The equations to solver maps to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    CALL Enters("SolverMapping_EquationsToSolverMapsInitialise",err,error,*999)

    IF(ASSOCIATED(equationsToSolverMaps)) THEN
      equationsToSolverMaps%equationsMatrixType=0
      equationsToSolverMaps%equationsMatrixNumber=0
      equationsToSolverMaps%solverMatrixNumber=0
      NULLIFY(equationsToSolverMaps%equationsMatrix)
      NULLIFY(equationsToSolverMaps%solverMatrix)
    ELSE
      CALL FlagError("Equations to solver map is not associated.",err,error,*999)
    ENDIF
        
    CALL Exits("SolverMapping_EquationsToSolverMapsInitialise")
    RETURN
999 CALL Errors("SolverMapping_EquationsToSolverMapsInitialise",err,error)    
    CALL Exits("SolverMapping_EquationsToSolverMapsInitialise")
    RETURN 1
   
  END SUBROUTINE SolverMapping_EquationsToSolverMapsInitialise

  !
  !================================================================================================================================
  !

  !>Finalises a equations set to solver matrix map and deallocates all memory.
  SUBROUTINE SolverMapping_EquationsToSolverInterfaceFinalise(equationsToSolverInterfaceMap,err,error,*)

    !Argument variables
    TYPE(EquationsToSolverMatrixMapsInterfaceType), INTENT(INOUT) :: equationsToSolverInterfaceMap !<The equations set to solver map interface to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL Enters("SolverMapping_EquationsToSolverInterfaceFinalise",err,error,*999)

    equationsToSolverInterfaceMap%interfaceConditionIndex=0
    NULLIFY(equationsToSolverInterfaceMap%interfaceCondition)
    equationsToSolverInterfaceMap%interfaceMatrixNumber=0
        
    CALL Exits("SolverMapping_EquationsToSolverInterfaceFinalise")
    RETURN
999 CALL Errors("SolverMapping_EquationsToSolverInterfaceFinalise",err,error)    
    CALL Exits("SolverMapping_EquationsToSolverInterfaceFinalise")
    RETURN 1
   
  END SUBROUTINE SolverMapping_EquationsToSolverInterfaceFinalise

  !
  !================================================================================================================================
  !

  !>Initialises a equations set to solver matrix interface map.
  SUBROUTINE SolverMapping_EquationsToSolverInterfaceInitialise(equationsToSolverInterfaceMap,err,error,*)

    !Argument variables
    TYPE(EquationsToSolverMatrixMapsInterfaceType), INTENT(OUT) :: equationsToSolverInterfaceMap !<The equations set to solver map interface to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL Enters("SolverMapping_EquationsToSolverInterfaceInitialise",err,error,*999)

    equationsToSolverInterfaceMap%interfaceConditionIndex=0
    NULLIFY(equationsToSolverInterfaceMap%interfaceCondition)
    equationsToSolverInterfaceMap%interfaceMatrixNumber=0
        
    CALL Exits("SolverMapping_EquationsToSolverInterfaceInitialise")
    RETURN
999 CALL Errors("SolverMapping_EquationsToSolverInterfaceInitialise",err,error)    
    CALL Exits("SolverMapping_EquationsToSolverInterfaceInitialise")
    RETURN 1
   
  END SUBROUTINE SolverMapping_EquationsToSolverInterfaceInitialise

  !
  !================================================================================================================================
  !

  !>Finalises a equations set to solver matrix map em and deallocates all memory.
  SUBROUTINE SolverMapping_EquationsToSolverMatrixMapsEmFinalise(equationsToSolverMatrixMapsEm,err,error,*)

    !Argument variables
    TYPE(EquationsToSolverMatrixMapsEmType), INTENT(INOUT) :: equationsToSolverMatrixMapsEm !<The equations set to solver matrix maps em to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: matrixIdx
    
    CALL Enters("SolverMapping_EquationsToSolverMatrixMapsEmFinalise",err,error,*999)
    
    IF(ALLOCATED(equationsToSolverMatrixMapsEm%equationsToSolverMatrixMaps)) THEN
      DO matrixIdx=1,SIZE(equationsToSolverMatrixMapsEm%equationsToSolverMatrixMaps,1)
        CALL SolverMapping_EquationsToSolverMapsFinalise(equationsToSolverMatrixMapsEm% &
          & equationsToSolverMatrixMaps(matrixIdx)%ptr,err,error,*999)        
      ENDDO !variableIdx
      DEALLOCATE(equationsToSolverMatrixMapsEm%equationsToSolverMatrixMaps)
    ENDIF
    
    CALL Exits("SolverMapping_EquationsToSolverMatrixMapsEmFinalise")
    RETURN
999 CALL Errors("SolverMapping_EquationsToSolverMatrixMapsEmFinalise",err,error)    
    CALL Exits("SolverMapping_EquationsToSolverMatrixMapsEmFinalise")
    RETURN 1
   
  END SUBROUTINE SolverMapping_EquationsToSolverMatrixMapsEmFinalise

  !
  !================================================================================================================================
  !

  !>Initialises an equations to solver matrix maps em.
  SUBROUTINE SolverMapping_EquationsToSolverMatrixMapsEmInitialise(equationsToSolverMatrixMapsEm,err,error,*)

    !Argument variables
    TYPE(EquationsToSolverMatrixMapsEmType), INTENT(OUT) :: equationsToSolverMatrixMapsEm !<The equations to solver matrix maps em to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    CALL Enters("SolverMapping_EquationsToSolverMatrixMapsEmInitialise",err,error,*999)

    equationsToSolverMatrixMapsEm%equationsMatrixNumber=0
    equationsToSolverMatrixMapsEm%numberOfSolverMatrices=0
        
    CALL Exits("SolverMapping_EquationsToSolverMatrixMapsEmInitialise")
    RETURN
999 CALL Errors("SolverMapping_EquationsToSolverMatrixMapsEmInitialise",err,error)    
    CALL Exits("SolverMapping_EquationsToSolverMatrixMapsEmInitialise")
    RETURN 1
   
  END SUBROUTINE SolverMapping_EquationsToSolverMatrixMapsEmInitialise

  !
  !================================================================================================================================
  !

  !>Finalises a equations set to solver matrix map jm and deallocates all memory.
  SUBROUTINE SolverMapping_EquationsToSolverMatrixMapsJmFinalise(equationsToSolverMatrixMapsJm,err,error,*)

    !Argument variables
    TYPE(JacobianToSolverMapPtrType), ALLOCATABLE, INTENT(IN) :: equationsToSolverMatrixMapsJm(:) !<equationsToSolverMatrixMapsJm(jacobianMatrixIdx). The equations set to solver matrix maps jm to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: equationsMatrixIdx

    CALL Enters("SolverMapping_EquationsToSolverMatrixMapsJmFinalise",err,error,*999)

    IF(ALLOCATED(equationsToSolverMatrixMapsJm)) THEN
      DO equationsMatrixIdx=1,SIZE(equationsToSolverMatrixMapsJm,1)
        CALL SolverMapping_JacobianToSolverMapFinalise(equationsToSolverMatrixMapsJm(equationsMatrixIdx)%ptr, &
          & err,error,*999)
      ENDDO !equationsMatrixIdx
    ENDIF

    CALL Exits("SolverMapping_EquationsToSolverMatrixMapsJmFinalise")
    RETURN
999 CALL Errors("SolverMapping_EquationsToSolverMatrixMapsJmFinalise",err,error)
    CALL Exits("SolverMapping_EquationsToSolverMatrixMapsJmFinalise")
    RETURN 1

  END SUBROUTINE SolverMapping_EquationsToSolverMatrixMapsJmFinalise

  !
  !================================================================================================================================
  !

  !>Finalises a equations set to solver matrix map sm and deallocates all memory.
  SUBROUTINE SolverMapping_EquationsToSolverMatrixMapsSmFinalise(equationsToSolverMatrixMapsSm,err,error,*)

    !Argument variables
    TYPE(EquationsToSolverMatrixMapsSmType), INTENT(INOUT) :: equationsToSolverMatrixMapsSm !<The equations set to solver matrix maps sm to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: matrixIdx,variableIdx
    
    CALL Enters("SolverMapping_EquationsToSolverMatrixMapsSmFinalise",err,error,*999)

    IF(ALLOCATED(equationsToSolverMatrixMapsSm%variableTypes)) DEALLOCATE(equationsToSolverMatrixMapsSm%variableTypes)
    IF(ALLOCATED(equationsToSolverMatrixMapsSm%variables)) DEALLOCATE(equationsToSolverMatrixMapsSm%variables)
    IF(ALLOCATED(equationsToSolverMatrixMapsSm%variableToSolverColMaps)) THEN
      DO variableIdx=1,SIZE(equationsToSolverMatrixMapsSm%variableToSolverColMaps,1)
        CALL SolverMapping_VariableToSolverColMapFinalise(equationsToSolverMatrixMapsSm% &
          & variableToSolverColMaps(variableIdx),err,error,*999)        
      ENDDO !variableIdx
      DEALLOCATE(equationsToSolverMatrixMapsSm%variableToSolverColMaps)
    ENDIF
    IF(ALLOCATED(equationsToSolverMatrixMapsSm%dynamicEquationsToSolverMatrixMaps)) THEN
      DO matrixIdx=1,SIZE(equationsToSolverMatrixMapsSm%dynamicEquationsToSolverMatrixMaps,1)
        CALL SolverMapping_EquationsToSolverMapsFinalise(equationsToSolverMatrixMapsSm% &
          & dynamicEquationsToSolverMatrixMaps(matrixIdx)%ptr,err,error,*999)        
      ENDDO !variableIdx
      DEALLOCATE(equationsToSolverMatrixMapsSm%dynamicEquationsToSolverMatrixMaps)
    ENDIF
    IF(ALLOCATED(equationsToSolverMatrixMapsSm%linearEquationsToSolverMatrixMaps)) THEN
      DO matrixIdx=1,SIZE(equationsToSolverMatrixMapsSm%linearEquationsToSolverMatrixMaps,1)
        CALL SolverMapping_EquationsToSolverMapsFinalise(equationsToSolverMatrixMapsSm% &
          & linearEquationsToSolverMatrixMaps(matrixIdx)%ptr,err,error,*999)        
      ENDDO !variableIdx
      DEALLOCATE(equationsToSolverMatrixMapsSm%linearEquationsToSolverMatrixMaps)
    ENDIF
    IF(ALLOCATED(equationsToSolverMatrixMapsSm%jacobianToSolverMatrixMaps)) THEN
      DO matrixIdx=1,SIZE(equationsToSolverMatrixMapsSm%jacobianToSolverMatrixMaps,1)
        CALL SolverMapping_JacobianToSolverMapFinalise(equationsToSolverMatrixMapsSm% &
          & jacobianToSolverMatrixMaps(matrixIdx)%ptr,err,error,*999)
      ENDDO
      DEALLOCATE(equationsToSolverMatrixMapsSm%jacobianToSolverMatrixMaps)
    ENDIF

    CALL Exits("SolverMapping_EquationsToSolverMatrixMapsSmFinalise")
    RETURN
999 CALL Errors("SolverMapping_EquationsToSolverMatrixMapsSmFinalise",err,error)    
    CALL Exits("SolverMapping_EquationsToSolverMatrixMapsSmFinalise")
    RETURN 1
   
  END SUBROUTINE SolverMapping_EquationsToSolverMatrixMapsSmFinalise

  !
  !================================================================================================================================
  !

  !>Initialises an equations to solver matrix maps sm.
  SUBROUTINE SolverMapping_EquationsToSolverMatrixMapsSmInitialise(equationsToSolverMatrixMapsSm,err,error,*)

    !Argument variables
    TYPE(EquationsToSolverMatrixMapsSmType), INTENT(OUT) :: equationsToSolverMatrixMapsSm !<The equations to solver matrix maps sm to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    CALL Enters("SolverMapping_EquationsToSolverMatrixMapsSmInitialise",err,error,*999)

    equationsToSolverMatrixMapsSm%solverMatrixNumber=0
    equationsToSolverMatrixMapsSm%numberOfVariables=0
    equationsToSolverMatrixMapsSm%numberOfDynamicEquationsMatrices=0
    equationsToSolverMatrixMapsSm%numberOfLinearEquationsMatrices=0
    equationsToSolverMatrixMapsSm%numberOfEquationsJacobians=0

    CALL Exits("SolverMapping_EquationsToSolverMatrixMapsSmInitialise")
    RETURN
999 CALL Errors("SolverMapping_EquationsToSolverMatrixMapsSmInitialise",err,error)    
    CALL Exits("SolverMapping_EquationsToSolverMatrixMapsSmInitialise")
    RETURN 1
   
  END SUBROUTINE SolverMapping_EquationsToSolverMatrixMapsSmInitialise

  !
  !================================================================================================================================
  !

  !>Finalises the solver mapping and deallocates all memory.
  SUBROUTINE SolverMapping_Finalise(solverMapping,err,error,*)

    !Argument variables
    TYPE(SolverMappingType), POINTER :: solverMapping !<A pointer to the solver mapping to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: equationsSetIdx,interfaceConditionIdx,rowIdx,solverMatrixIdx

    CALL Enters("SolverMapping_Finalise",err,error,*999)

    IF(ASSOCIATED(solverMapping)) THEN
      CALL SolverMapping_VariablesFinalise(solverMapping%RowVariablesList,err,error,*999)
      IF(ALLOCATED(solverMapping%columnVariablesList)) THEN
        DO solverMatrixIdx=1,SIZE(SolverMapping%columnVariablesList,1)
          CALL SolverMapping_VariablesFinalise(solverMapping%columnVariablesList(solverMatrixIdx),err,error,*999)
        ENDDO ! solverMatrixIdx
        DEALLOCATE(solverMapping%columnVariablesList)
      ENDIF
      IF(ALLOCATED(solverMapping%equationsSets)) DEALLOCATE(solverMapping%equationsSets)        
      IF(ALLOCATED(solverMapping%equationsSetToSolverMap)) THEN
        DO equationsSetIdx=1,SIZE(solverMapping%equationsSetToSolverMap,1)
          CALL SolverMapping_EquationsSetToSolverMapFinalise(solverMapping%equationsSetToSolverMap(equationsSetIdx), &
            & err,error,*999)
        ENDDO !equationsSetIdx
        DEALLOCATE(solverMapping%equationsSetToSolverMap)
      ENDIF
      IF(ALLOCATED(solverMapping%interfaceConditions)) DEALLOCATE(solverMapping%interfaceConditions)
      IF(ALLOCATED(solverMapping%interfaceConditionToSolverMap)) THEN
        DO interfaceConditionIdx=1,SIZE(solverMapping%interfaceConditionToSolverMap,1)
          CALL SolverMapping_InterfaceConditionToSolverMapFinalise(solverMapping%interfaceConditionToSolverMap( &
            & interfaceConditionIdx),err,error,*999)
        ENDDO !interfaceConditionIdx
        DEALLOCATE(solverMapping%interfaceConditionToSolverMap)
      ENDIF
      IF(ALLOCATED(solverMapping%solverColToEquationsColsMap)) THEN
        DO solverMatrixIdx=1,SIZE(solverMapping%solverColToEquationsColsMap,1)
          CALL SolverMapping_SolverColToEquationsMapsFinalise(solverMapping%solverColToEquationsColsMap(solverMatrixIdx), &
            & err,error,*999)
        ENDDO !solverMatrixIdx
        DEALLOCATE(solverMapping%solverColToEquationsColsMap)
      ENDIF
      IF(ALLOCATED(solverMapping%solverRowToEquationsRowsMap)) THEN
        DO rowIdx=1,SIZE(solverMapping%solverRowToEquationsRowsMap,1)
          CALL SolverMapping_SolverRowToEquationsMapsFinalise(solverMapping%solverRowToEquationsRowsMap(rowIdx),err,error,*999)
        ENDDO !rowIdx
        DEALLOCATE(solverMapping%solverRowToEquationsRowsMap)
      ENDIF
      CALL DOMAIN_MAPPINGS_MAPPING_FINALISE(solverMapping%rowDofsMapping,err,error,*999)
      CALL SolverMapping_CreateValuesCacheFinalise(solverMapping%createValuesCache,err,error,*999)
      DEALLOCATE(solverMapping)
    ENDIF
       
    CALL Exits("SolverMapping_Finalise")
    RETURN
999 CALL Errors("SolverMapping_Finalise",err,error)
    CALL Exits("SolverMapping_Finalise")
    RETURN 1
  END SUBROUTINE SolverMapping_Finalise

  !
  !================================================================================================================================
  !

  !>Initialises the solver mapping and deallocates all memory.
  SUBROUTINE SolverMapping_Initialise(solverEquations,err,error,*)

    !Argument variables
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: solverEquations !<A pointer to the solver equations to initialise the solver mapping on.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    CALL Enters("SolverMapping_Initialise",err,error,*998)

    IF(ASSOCIATED(solverEquations)) THEN
      IF(ASSOCIATED(solverEquations%SOLVER_MAPPING)) THEN
        CALL FlagError("Solver equations solver mapping is already associated.",err,error,*998)
      ELSE
        ALLOCATE(solverEquations%SOLVER_MAPPING,STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate solver equations solver mapping.",err,error,*999)
        solverEquations%SOLVER_MAPPING%solverEquations=>solverEquations
        solverEquations%SOLVER_MAPPING%solverMappingFinished=.FALSE.
        solverEquations%SOLVER_MAPPING%numberOfSolverMatrices=1
        solverEquations%SOLVER_MAPPING%numberOfRows=0
        solverEquations%SOLVER_MAPPING%numberOfGlobalRows=0
        solverEquations%SOLVER_MAPPING%numberOfEquationsSets=0
        solverEquations%SOLVER_MAPPING%numberOfInterfaceConditions=0
        NULLIFY(solverEquations%SOLVER_MAPPING%rowDofsMapping)
        NULLIFY(solverEquations%SOLVER_MAPPING%createValuesCache)
        CALL SolverMapping_CreateValuesCacheInitialise(solverEquations%SOLVER_MAPPING,err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Solver equations is not associated.",err,error,*998)
    ENDIF
    
    CALL Exits("SolverMapping_Initialise")
    RETURN
999 CALL SolverMapping_Finalise(solverEquations%SOLVER_MAPPING,dummyErr,dummyError,*998)
998 CALL Errors("SolverMapping_Initialise",err,error)
    CALL Exits("SolverMapping_Initialise")
    RETURN 1
  END SUBROUTINE SolverMapping_Initialise

  !
  !================================================================================================================================
  !

  !>Adds an interface condition to a solver mapping
  SUBROUTINE SolverMapping_InterfaceConditionAdd(solverMapping,interfaceCondition,interfaceConditionIndex,err,error,*)

    !Argument variables
    TYPE(SolverMappingType), POINTER :: solverMapping !<A pointer the solver mapping to add the interface condition to
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: interfaceCondition !<A pointer to the interface condition to add
    INTEGER(INTG), INTENT(OUT) :: interfaceConditionIndex !<On exit, the index of the interface condition in the solver mapping
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: equationsMatrixIdx,equationsSetIdx,interfaceConditionIdx,interfaceMatrixIdx,listItem(2), &
      & numberOfInterfaceMatrices
    LOGICAL :: equationsSetFound,variableFound
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: dependentVariable
    TYPE(INTERFACE_DEPENDENT_TYPE), POINTER :: interfaceDependent
    TYPE(INTERFACE_EQUATIONS_TYPE), POINTER :: interfaceEquations
    TYPE(INTERFACE_MAPPING_TYPE), POINTER :: interfaceMapping
    TYPE(INTERFACE_CONDITION_ptr_TYPE), ALLOCATABLE :: newInterfaceConditions(:)
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: solverEquations
    TYPE(VARYING_STRING) :: localError

    CALL Enters("SolverMapping_InterfaceConditionAdd",err,error,*999)

    InterfaceConditionIndex=0
    IF(ASSOCIATED(solverMapping)) THEN
      IF(solverMapping%solverMappingFinished) THEN
        CALL FlagError("Solver mapping has been finished.",err,error,*999)
      ELSE
        solverEquations=>solverMapping%solverEquations
        IF(ASSOCIATED(solverEquations)) THEN
          IF(ASSOCIATED(interfaceCondition)) THEN
            IF(interfaceCondition%INTERFACE_CONDITION_FINISHED) THEN
              interfaceEquations=>interfaceCondition%INTERFACE_EQUATIONS
              IF(ASSOCIATED(interfaceEquations)) THEN
                interfaceMapping=>interfaceEquations%INTERFACE_MAPPING
                IF(ASSOCIATED(interfaceMapping)) THEN                
                  IF(ASSOCIATED(solverMapping%createValuesCache)) THEN
                    !Check that the interface variables are already part of an added equations set.
                    SELECT CASE(interfaceCondition%method)
                    CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD,INTERFACE_CONDITION_PENALTY_METHOD)
                      interfaceDependent=>interfaceCondition%DEPENDENT
                      IF(ASSOCIATED(interfaceDependent)) THEN
                        SELECT CASE(interfaceCondition%method)
                        CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD)
                          numberOfInterfaceMatrices=interfaceMapping%NUMBER_OF_INTERFACE_MATRICES
                        CASE(INTERFACE_CONDITION_PENALTY_METHOD)
                          numberOfInterfaceMatrices=interfaceMapping%NUMBER_OF_INTERFACE_MATRICES-1
                        CASE DEFAULT
                          localError="The interface condition method of "// &
                            & TRIM(NumberToVstring(interfaceCondition%method,"*",err,error))//" is invalid."
                          CALL FlagError(localError,err,error,*999)
                        ENDSELECT
                        DO interfaceMatrixIdx=1,numberOfInterfaceMatrices
                          equationsSet=>interfaceMapping%INTERFACE_MATRIX_ROWS_TO_VAR_MAPS(interfaceMatrixIdx)%EQUATIONS_SET
                          IF(ASSOCIATED(equationsSet)) THEN
                            equationsSetFound=.FALSE.
                            DO equationsSetIdx=1,solverMapping%numberOfEquationsSets
                              IF(ASSOCIATED(equationsSet,solverMapping%equationsSets(equationsSetIdx)%ptr)) THEN
                                equationsSetFound=.TRUE.
                                EXIT
                              ENDIF
                            ENDDO !equationsSetIdx
                            IF(equationsSetFound) THEN
                              !See if the variable is in the equations set.
                              dependentVariable=>interfaceMapping%INTERFACE_MATRIX_ROWS_TO_VAR_MAPS(interfaceMatrixIdx)%variable
                              IF(ASSOCIATED(dependentVariable)) THEN
                                variableFound=.FALSE.
                                !Check dynamic variables
                                IF(solverMapping%createValuesCache%dynamicVariableType(equationsSetIdx)== &
                                  & dependentVariable%VARIABLE_TYPE) THEN
                                  variableFound=.TRUE.
                                ELSE
                                  !Check linear matrices. Just check for solver matrix 1 and the moment
                                  DO equationsMatrixIdx=1,solverMapping%createValuesCache% &
                                    & matrixVariableTypes(0,equationsSetIdx,1)
                                    IF(solverMapping%createValuesCache% &
                                      & matrixVariableTypes(equationsMatrixIdx,equationsSetIdx,1)== &
                                      & dependentVariable%VARIABLE_TYPE) THEN
                                      variableFound=.TRUE.
                                      EXIT
                                    ENDIF
                                  ENDDO !equationsMatrixIdx
                                  IF(.NOT.variableFound) THEN
                                    !Check residual variable type
                                    DO equationsMatrixIdx=1,solverMapping%createValuesCache%residualVariableTypes(0, &
                                      & equationsSetIdx)
                                      IF(solverMapping%createValuesCache%residualVariableTypes(equationsMatrixIdx, &
                                        & equationsSetIdx)==dependentVariable%VARIABLE_TYPE) THEN
                                        variableFound=.TRUE.
                                        EXIT
                                      ENDIF
                                    ENDDO
                                  ENDIF
                                ENDIF
                                IF(variableFound) THEN
                                  !Add in interface condition to equations set (just for solver matrix 1 at the moment)
                                  listItem(1)=solverMapping%numberOfInterfaceConditions+1
                                  listItem(2)=interfaceMatrixIdx
                                  CALL List_ItemAdd(solverMapping%createValuesCache%interfaceIndices(equationsSetIdx)%ptr, &
                                    & listItem,err,error,*999)                                  
                                ELSE
                                  localError="The dependent variable associated with interface matrix number "// &
                                    & TRIM(NumberToVString(interfaceMatrixIdx,"*",err,error))// &
                                    & " is not mapped to the solver equations."
                                  CALL FlagError(localError,err,error,*999)
                                ENDIF
                              ELSE
                                localError="The dependent variable associated with interface matrix number "// &
                                  & TRIM(NumberToVString(interfaceMatrixIdx,"*",err,error))//" is not associated."
                                CALL FlagError(localError,err,error,*999)
                              ENDIF
                            ELSE
                              localError="The equations set for the dependent variable associated with interface "// &
                                & "matrix number "//TRIM(NumberToVString(interfaceMatrixIdx,"*",err,error))// &
                                & " has not been added to the solver equations."
                              CALL FlagError(localError,err,error,*999)
                            ENDIF
                          ELSE
                            localError="Equations set is not associated for interface matrix number "// &
                              & TRIM(NumberToVString(interfaceMatrixIdx,"*",err,error))//"."
                            CALL FlagError(localError,err,error,*999)
                          ENDIF
                        ENDDO !interfaceMatrixIdx
                      ELSE
                        CALL FlagError("Interface condition dependent is not associated.",err,error,*999)
                      ENDIF
                    CASE(INTERFACE_CONDITION_AUGMENTED_LAGRANGE_METHOD)
                      CALL FlagError("Not implemented.",err,error,*999)
                    CASE(INTERFACE_CONDITION_POINT_TO_POINT_METHOD)
                      CALL FlagError("Not implemented.",err,error,*999)
                    CASE DEFAULT
                      localError="The interface condition method of "// &
                        & TRIM(NumberToVString(interfaceCondition%method,"*",err,error))//" is invalid."
                      CALL FlagError(localError,err,error,*999)
                    END SELECT
                    IF(solverMapping%numberOfInterfaceConditions>0) THEN
                      ALLOCATE(newInterfaceConditions(solverMapping%numberOfInterfaceConditions+1),STAT=err)
                      IF(err/=0) CALL FlagError("Could not allocate new interface conditions.",err,error,*999)                   
                      DO interfaceConditionIdx=1,solverMapping%numberOfInterfaceConditions
                        newInterfaceConditions(interfaceConditionIdx)%ptr=> &
                          & solverMapping%interfaceConditions(interfaceConditionIdx)%ptr
                      ENDDO !interfaceConditionIdx
                      newInterfaceConditions(solverMapping%numberOfInterfaceConditions+1)%ptr=>interfaceCondition
                    ELSE
                      ALLOCATE(newInterfaceConditions(1),STAT=err)
                      IF(err/=0) CALL FlagError("Could not allocate new interface conditions.",err,error,*999)
                      newInterfaceConditions(1)%ptr=>interfaceCondition
                    ENDIF
                    CALL MOVE_ALLOC(newInterfaceConditions,solverMapping%interfaceConditions)
                    solverMapping%numberOfInterfaceConditions=solverMapping%numberOfInterfaceConditions+1
                    interfaceConditionIndex=solverMapping%numberOfInterfaceConditions

!!TODO: SORT OUT LAGRANGE FIELD VARIABLE
                    CALL SolverMapping_CreateValuesCacheInterfaceVarListAdd(solverMapping,solverMapping%createValuesCache% &
                      & interfaceColVariablesList(1)%ptr,interfaceConditionIndex,1,err,error,*999)
                    
                  ELSE
                    CALL FlagError("Solvers mapping create values cache is not associated.",err,error,*999)
                  ENDIF
                ELSE
                  CALL FlagError("Interface equations mapping is not associated.",err,error,*999)
                ENDIF
              ELSE
                CALL FlagError("Interface condition interface equations is not associated.",err,error,*999)
              ENDIF
            ELSE
              CALL FlagError("Interface condition has not been finished.",err,error,*999)
            ENDIF
          ELSE
            CALL FlagError("Interface condition is not associated.",err,error,*999)
          ENDIF
        ELSE
          CALL FlagError("Solver mapping solver is not associated.",err,error,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Solver mapping is not associated.",err,error,*999)
    ENDIF
        
    CALL Exits("SolverMapping_InterfaceConditionAdd")
    RETURN
999 IF(ALLOCATED(newInterfaceConditions)) DEALLOCATE(newInterfaceConditions)
    CALL Errors("SolverMapping_InterfaceConditionAdd",err,error)    
    CALL Exits("SolverMapping_InterfaceConditionAdd")
    RETURN 1
   
  END SUBROUTINE SolverMapping_InterfaceConditionAdd

  !
  !================================================================================================================================
  !

  !>Finalises an interface condition to solver map and deallocates all memory.
  SUBROUTINE SolverMapping_InterfaceConditionToSolverMapFinalise(interfaceConditionToSolverMap,err,error,*)

    !Argument variables
    TYPE(interfaceConditionToSolverMapType), INTENT(INOUT) :: interfaceConditionToSolverMap !<The interface condition to solver map to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: columnIdx,equationsSetIdx,interfaceMatrixIdx,solverMatrixIdx
    
    CALL Enters("SolverMapping_InterfaceConditionToSolverMapFinalise",err,error,*999)

    interfaceConditionToSolverMap%interfaceConditionIndex=0
    NULLIFY(interfaceConditionToSolverMap%solverMapping)
    NULLIFY(interfaceConditionToSolverMap%interfaceEquations)
    interfaceConditionToSolverMap%numberOfEquationsSets=0
    IF(ALLOCATED(interfaceConditionToSolverMap%interfaceToSolverMatrixMapsEquations)) THEN
      DO equationsSetIdx=1,SIZE(interfaceConditionToSolverMap%interfaceToSolverMatrixMapsEquations,1)
        CALL SolverMapping_InterfaceToSolverEquationsFinalise(interfaceConditionToSolverMap% &
          & interfaceToSolverMatrixMapsEquations(equationsSetIdx),err,error,*999)
      ENDDO !equationsSetIdx
      DEALLOCATE(interfaceConditionToSolverMap%interfaceToSolverMatrixMapsEquations)
    ENDIF
    IF(ALLOCATED(interfaceConditionToSolverMap%interfaceToSolverMatrixMapsSm)) THEN
      DO solverMatrixIdx=1,SIZE(interfaceConditionToSolverMap%interfaceToSolverMatrixMapsSm,1)
        CALL SolverMapping_InterfaceToSolverMatrixMapsSmFinalise(interfaceConditionToSolverMap% &
          & interfaceToSolverMatrixMapsSm(solverMatrixIdx),err,error,*999)
      ENDDO !solverMatrixIdx
      DEALLOCATE(interfaceConditionToSolverMap%interfaceToSolverMatrixMapsSm)
    ENDIF
    IF(ALLOCATED(interfaceConditionToSolverMap%interfaceToSolverMatrixMapsIm)) THEN
      DO interfaceMatrixIdx=1,SIZE(interfaceConditionToSolverMap%interfaceToSolverMatrixMapsIm,1)
        CALL SolverMapping_InterfaceToSolverMatrixMapsImFinalise(interfaceConditionToSolverMap% &
          & interfaceToSolverMatrixMapsIm(interfaceMatrixIdx),err,error,*999)
      ENDDO !interfaceMatrixIdx
      DEALLOCATE(interfaceConditionToSolverMap%interfaceToSolverMatrixMapsIm)
    ENDIF
    IF(ALLOCATED(interfaceConditionToSolverMap%interfaceColumnToSolverRowsMaps)) THEN
      DO columnIdx=1,SIZE(interfaceConditionToSolverMap%interfaceColumnToSolverRowsMaps,1)
        CALL SolverMapping_InterfaceColToSolverRowsMapFinalise(interfaceConditionToSolverMap% &
          & interfaceColumnToSolverRowsMaps(columnIdx),err,error,*999)
      ENDDO !columnIdx
      DEALLOCATE(interfaceConditionToSolverMap%interfaceColumnToSolverRowsMaps)
    ENDIF
        
    CALL Exits("SolverMapping_InterfaceConditionToSolverMapFinalise")
    RETURN
999 CALL Errors("SolverMapping_InterfaceConditionToSolverMapFinalise",err,error)    
    CALL Exits("SolverMapping_InterfaceConditionToSolverMapFinalise")
    RETURN 1
   
  END SUBROUTINE SolverMapping_InterfaceConditionToSolverMapFinalise

  !
  !================================================================================================================================
  !

  !>Initialises an interface condition to solver map.
  SUBROUTINE SolverMapping_InterfaceConditionToSolverMapInitialise(interfaceConditionToSolverMap,err,error,*)

    !Argument variables
    TYPE(interfaceConditionToSolverMapType), INTENT(OUT) :: interfaceConditionToSolverMap !<The interface condition to solver map to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL Enters("SolverMapping_InterfaceConditionToSolverMapInitialise",err,error,*999)

    interfaceConditionToSolverMap%interfaceConditionIndex=0
    NULLIFY(interfaceConditionToSolverMap%solverMapping)
    NULLIFY(interfaceConditionToSolverMap%interfaceEquations)
    interfaceConditionToSolverMap%numberOfEquationsSets=0
    
    CALL Exits("SolverMapping_InterfaceConditionToSolverMapInitialise")
    RETURN
999 CALL Errors("SolverMapping_InterfaceConditionToSolverMapInitialise",err,error)    
    CALL Exits("SolverMapping_InterfaceConditionToSolverMapInitialise")
    RETURN 1
   
  END SUBROUTINE SolverMapping_InterfaceConditionToSolverMapInitialise

  !
  !================================================================================================================================
  !

  !>Finalises an interface to solver matrix maps and deallocates all memory.
  SUBROUTINE SolverMapping_InterfaceToSolverMapsFinalise(interfaceToSolverMaps,err,error,*)

    !Argument variables
    TYPE(InterfaceToSolverMapsType), POINTER :: interfaceToSolverMaps !<The interface to solver maps to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: rowIdx
    
    CALL Enters("SolverMapping_InterfaceToSolverMapsFinalise",err,error,*999)

    IF(ASSOCIATED(interfaceToSolverMaps)) THEN
      IF(ALLOCATED(interfaceToSolverMaps%interfaceRowToSolverColsMap)) THEN
        DO rowIdx=1,SIZE(interfaceToSolverMaps%interfaceRowToSolverColsMap,1)
          CALL SolverMapping_EquationsColToSolverColsMapFinalise(interfaceToSolverMaps% &
            & interfaceRowToSolverColsMap(rowIdx),err,error,*999)
        ENDDO !rowIdx
        DEALLOCATE(interfaceToSolverMaps%interfaceRowToSolverColsMap)
      ENDIF
    ENDIF
        
    CALL Exits("SolverMapping_InterfaceToSolverMapsFinalise")
    RETURN
999 CALL Errors("SolverMapping_InterfaceToSolverMapsFinalise",err,error)    
    CALL Exits("SolverMapping_InterfaceToSolverMapsFinalise")
    RETURN 1
   
  END SUBROUTINE SolverMapping_InterfaceToSolverMapsFinalise

  !
  !================================================================================================================================
  !

  !>Initialises an interface to solver maps
  SUBROUTINE SolverMapping_InterfaceConditionToSolverMapsInitialise(interfaceToSolverMaps,err,error,*)

    !Argument variables
    TYPE(InterfaceToSolverMapsType), POINTER :: interfaceToSolverMaps !<The interface to solver maps to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    CALL Enters("SolverMapping_InterfaceConditionToSolverMapsInitialise",err,error,*999)

    IF(ASSOCIATED(interfaceToSolverMaps)) THEN
      interfaceToSolverMaps%interfaceMatrixNumber=0
      interfaceToSolverMaps%solverMatrixNumber=0
      NULLIFY(interfaceToSolverMaps%interfaceMatrix)
      NULLIFY(interfaceToSolverMaps%solverMatrix)
    ELSE
      CALL FlagError("Interface to solver maps is not associated.",err,error,*999)
    ENDIF
        
    CALL Exits("SolverMapping_InterfaceConditionToSolverMapsInitialise")
    RETURN
999 CALL Errors("SolverMapping_InterfaceConditionToSolverMapsInitialise",err,error)    
    CALL Exits("SolverMapping_InterfaceConditionToSolverMapsInitialise")
    RETURN 1
   
  END SUBROUTINE SolverMapping_InterfaceConditionToSolverMapsInitialise

  !
  !================================================================================================================================
  !

  !>Finalises an interface to solver matrix equations map and deallocates all memory.
  SUBROUTINE SolverMapping_InterfaceToSolverEquationsFinalise(interfaceToSolverEquationsMaps,err,error,*)

    !Argument variables
    TYPE(interfaceToSolverMatrixMapsEquationsType), INTENT(INOUT) :: interfaceToSolverEquationsMaps !<The interface to solver equations map to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    CALL Enters("SolverMapping_InterfaceToSolverEquationsFinalise",err,error,*999)

    interfaceToSolverEquationsMaps%equationsSetIndex=0
    NULLIFY(interfaceToSolverEquationsMaps%equationsSet)
    interfaceToSolverEquationsMaps%interfaceMatrixIndex=0
         
    CALL Exits("SolverMapping_InterfaceToSolverEquationsFinalise")
    RETURN
999 CALL Errors("SolverMapping_InterfaceToSolverEquationsFinalise",err,error)    
    CALL Exits("SolverMapping_InterfaceToSolverEquationsFinalise")
    RETURN 1
   
  END SUBROUTINE SolverMapping_InterfaceToSolverEquationsFinalise

  !
  !================================================================================================================================
  !

  !>Initialises an interface to solver matrix equations map.
  SUBROUTINE SolverMapping_InterfaceToSolverEquationsInitialise(interfacetoSolverEquationsMaps,err,error,*)

    !Argument variables
    TYPE(interfaceToSolverMatrixMapsEquationsType), INTENT(OUT) :: interfacetoSolverEquationsMaps !<The interface to solver equations map to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    CALL Enters("SolverMapping_InterfaceToSolverEquationsInitialise",err,error,*999)

    interfacetoSolverEquationsMaps%equationsSetIndex=0
    NULLIFY(interfacetoSolverEquationsMaps%equationsSet)
    interfacetoSolverEquationsMaps%interfaceMatrixIndex=0
         
    CALL Exits("SolverMapping_InterfaceToSolverEquationsInitialise")
    RETURN
999 CALL Errors("SolverMapping_InterfaceToSolverEquationsInitialise",err,error)    
    CALL Exits("SolverMapping_InterfaceToSolverEquationsInitialise")
    RETURN 1
   
  END SUBROUTINE SolverMapping_InterfaceToSolverEquationsInitialise

  !
  !================================================================================================================================
  !

  !>Finalises a interface column to solver row map and deallocates all memory.
  SUBROUTINE SolverMapping_InterfaceColToSolverRowsMapFinalise(interfaceColumnToSolverRowsMap,err,error,*)

    !Argument variables
    TYPE(InterfaceColumnToSolverRowsMapType), INTENT(INOUT) :: interfaceColumnToSolverRowsMap !<The interface column to solver row map to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL Enters("SolverMapping_InterfaceColToSolverRowsMapFinalise",err,error,*999)

    interfaceColumnToSolverRowsMap%numberOfSolverRows=0
    interfaceColumnToSolverRowsMap%solverRow=0
    interfaceColumnToSolverRowsMap%couplingCoefficient=0.0_DP
         
    CALL Exits("SolverMapping_InterfaceColToSolverRowsMapFinalise")
    RETURN
999 CALL Errors("SolverMapping_InterfaceColToSolverRowsMapFinalise",err,error)    
    CALL Exits("SolverMapping_InterfaceColToSolverRowsMapFinalise")
    RETURN 1
   
  END SUBROUTINE SolverMapping_InterfaceColToSolverRowsMapFinalise

  !
  !================================================================================================================================
  !

  !>Initialises am interface column to solver row map.
  SUBROUTINE SolverMapping_InterfaceColToSolverRowsMapInitialise(interfaceColumnToSolverRowsMap,err,error,*)

    !Argument variables
    TYPE(InterfaceColumnToSolverRowsMapType), INTENT(OUT) :: interfaceColumnToSolverRowsMap !<The interface column to solver row map to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    CALL Enters("SolverMapping_InterfaceColToSolverRowsMapInitialise",err,error,*999)

    interfaceColumnToSolverRowsMap%numberOfSolverRows=0
    interfaceColumnToSolverRowsMap%solverRow=0
    interfaceColumnToSolverRowsMap%couplingCoefficient=0.0_DP
        
    CALL Exits("SolverMapping_InterfaceColToSolverRowsMapInitialise")
    RETURN
999 CALL Errors("SolverMapping_InterfaceColToSolverRowsMapInitialise",err,error)    
    CALL Exits("SolverMapping_InterfaceColToSolverRowsMapInitialise")
    RETURN 1
   
  END SUBROUTINE SolverMapping_InterfaceColToSolverRowsMapInitialise

  !
  !================================================================================================================================
  !

  !>Finalises a interface row to solver row map and deallocates all memory.
  SUBROUTINE SolverMapping_InterfaceRowToSolverRowsMapFinalise(interfaceRowToSolverRowsMap,err,error,*)

    !Argument variables
    TYPE(InterfaceRowToSolverRowsMapType), INTENT(OUT) :: interfaceRowToSolverRowsMap !<The interface row to solver row map to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    CALL Enters("SolverMapping_InterfaceRowToSolverRowsMapFinalise",err,error,*999)

    interfaceRowToSolverRowsMap%numberOfSolverRows=0
    interfaceRowToSolverRowsMap%solverRow=0
    interfaceRowToSolverRowsMap%couplingCoefficient=0.0_DP
        
    CALL Exits("SolverMapping_InterfaceRowToSolverRowsMapFinalise")
    RETURN
999 CALL Errors("SolverMapping_InterfaceRowToSolverRowsMapFinalise",err,error)    
    CALL Exits("SolverMapping_InterfaceRowToSolverRowsMapFinalise")
    RETURN 1
   
  END SUBROUTINE SolverMapping_InterfaceRowToSolverRowsMapFinalise

  !
  !================================================================================================================================
  !

  !>Initialises am interface row to solver row map.
  SUBROUTINE SolverMapping_InterfaceRowToSolverRowsMapInitialise(interfaceRowToSolverRowsMap,err,error,*)

    !Argument variables
    TYPE(InterfaceRowToSolverRowsMapType), INTENT(OUT) :: interfaceRowToSolverRowsMap !<The interface row to solver row map to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    CALL Enters("SolverMapping_InterfaceRowToSolverRowsMapInitialise",err,error,*999)

    interfaceRowToSolverRowsMap%numberOfSolverRows=0
    interfaceRowToSolverRowsMap%solverRow=0
    interfaceRowToSolverRowsMap%couplingCoefficient=0.0_DP
        
    CALL Exits("SolverMapping_InterfaceRowToSolverRowsMapInitialise")
    RETURN
999 CALL Errors("SolverMapping_InterfaceRowToSolverRowsMapInitialise",err,error)    
    CALL Exits("SolverMapping_InterfaceRowToSolverRowsMapInitialise")
    RETURN 1
   
  END SUBROUTINE SolverMapping_InterfaceRowToSolverRowsMapInitialise

  !
  !================================================================================================================================
  !

  !>Finalises an interface to solver matrix map im and deallocates all memory.
  SUBROUTINE SolverMapping_InterfaceToSolverMatrixMapsImFinalise(interfaceToSolverMatrixMapsIm,err,error,*)

    !Argument variables
    TYPE(InterfaceToSolverMatrixMapsImType), INTENT(INOUT) :: interfaceToSolverMatrixMapsIm !<The interface to solver matrix maps Im to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: rowIdx,solverMatrixIdx
    
    CALL Enters("SolverMapping_InterfaceToSolverMatrixMapsImFinalise",err,error,*999)

    interfaceToSolverMatrixMapsIm%interfaceMatrixNumber=0
    interfaceToSolverMatrixMapsIm%numberOfSolverMatrices=0
    IF(ALLOCATED(interfaceToSolverMatrixMapsIm%interfaceToSolverMatrixMaps)) THEN
      DO solverMatrixIdx=1,SIZE(interfaceToSolverMatrixMapsIm%interfaceToSolverMatrixMaps,1)
        CALL SolverMapping_InterfaceToSolverMapsFinalise(interfaceToSolverMatrixMapsIm% &
          & interfaceToSolverMatrixMaps(solverMatrixIdx)%ptr,err,error,*999)        
      ENDDO !solverMatrixIdx
      DEALLOCATE(interfaceToSolverMatrixMapsIm%interfaceToSolverMatrixMaps)
    ENDIF
    IF(ALLOCATED(interfaceToSolverMatrixMapsIm%interfaceRowToSolverRowsMap)) THEN
      DO rowIdx=1,SIZE(interfaceToSolverMatrixMapsIm%interfaceRowToSolverRowsMap,1)
        CALL SolverMapping_InterfaceRowToSolverRowsMapFinalise(interfaceToSolverMatrixMapsIm% &
          interfaceRowToSolverRowsMap(rowIdx),err,error,*999)
      ENDDO !rowIdx
      DEALLOCATE(interfaceToSolverMatrixMapsIm%interfaceRowToSolverRowsMap)
    ENDIF
    
    CALL Exits("SolverMapping_InterfaceToSolverMatrixMapsImFinalise")
    RETURN
999 CALL Errors("SolverMapping_InterfaceToSolverMatrixMapsImFinalise",err,error)    
    CALL Exits("SolverMapping_InterfaceToSolverMatrixMapsImFinalise")
    RETURN 1
   
  END SUBROUTINE SolverMapping_InterfaceToSolverMatrixMapsImFinalise

  !
  !================================================================================================================================
  !

  !>Initialises an interface to solver matrix maps im.
  SUBROUTINE SolverMapping_InterfaceToSolverMatrixMapsImInitialise(interfaceToSolverMatrixMapsIm,err,error,*)

    !Argument variables
    TYPE(InterfaceToSolverMatrixMapsImType), INTENT(OUT) :: interfaceToSolverMatrixMapsIm !<The interface to solver matrix maps Im to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    CALL Enters("SolverMapping_InterfaceToSolverMatrixMapsImInitialise",err,error,*999)

    interfaceToSolverMatrixMapsIm%interfaceMatrixNumber=0
    interfaceToSolverMatrixMapsIm%numberOfSolverMatrices=0
    
    CALL Exits("SolverMapping_InterfaceToSolverMatrixMapsImInitialise")
    RETURN
999 CALL Errors("SolverMapping_InterfaceToSolverMatrixMapsImInitialise",err,error)    
    CALL Exits("SolverMapping_InterfaceToSolverMatrixMapsImInitialise")
    RETURN 1
   
  END SUBROUTINE SolverMapping_InterfaceToSolverMatrixMapsImInitialise

  !
  !================================================================================================================================
  !

  !>Finalises an interface to solver matrix map sm and deallocates all memory.
  SUBROUTINE SolverMapping_InterfaceToSolverMatrixMapsSmFinalise(interfaceToSolverMatrixMapsSm,err,error,*)

    !Argument variables
    TYPE(InterfaceToSolverMatrixMapsSmType), INTENT(INOUT) :: interfaceToSolverMatrixMapsSm !<The interface to solver matrix maps sm to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: columnIdx,interfaceMatrixIdx
    
    CALL Enters("SolverMapping_InterfaceToSolverMatrixMapsSmFinalise",err,error,*999)

    interfaceToSolverMatrixMapsSm%solverMatrixNumber=0
    interfaceToSolverMatrixMapsSm%lagrangeVariableType=0
    NULLIFY(interfaceToSolverMatrixMapsSm%lagrangeVariable)
    CALL SolverMapping_VariableToSolverColMapFinalise(interfaceToSolverMatrixMapsSm% &
      & lagrangeVariableToSolverColMap,err,error,*999)
    interfaceToSolverMatrixMapsSm%numberOfDependentVariables=0
    IF(ALLOCATED(interfaceToSolverMatrixMapsSm%dependentVariableTypes)) &
      & DEALLOCATE(interfaceToSolverMatrixMapsSm%dependentVariableTypes)
    IF(ALLOCATED(interfaceToSolverMatrixMapsSm%dependentVariables)) &
      & DEALLOCATE(interfaceToSolverMatrixMapsSm%dependentVariables)
    IF(ALLOCATED(interfaceToSolverMatrixMapsSm%dependentVariableToSolverColMaps)) THEN
      DO interfaceMatrixIdx=1,SIZE(interfaceToSolverMatrixMapsSm%dependentVariableToSolverColMaps,1)
        CALL SolverMapping_VariableToSolverColMapFinalise(interfaceToSolverMatrixMapsSm% &
          & dependentVariableToSolverColMaps(interfaceMatrixIdx),err,error,*999)
      ENDDO !interfaceMatrixIdx
      DEALLOCATE(interfaceToSolverMatrixMapsSm%dependentVariableToSolverColMaps)
    ENDIF
    interfaceToSolverMatrixMapsSm%numberOfInterfaceMatrices=0
    IF(ALLOCATED(interfaceToSolverMatrixMapsSm%interfaceEquationsToSolverMatrixMaps)) THEN
      DO interfaceMatrixIdx=1,SIZE(interfaceToSolverMatrixMapsSm%interfaceEquationsToSolverMatrixMaps,1)
        CALL SolverMapping_InterfaceToSolverMapsFinalise(interfaceToSolverMatrixMapsSm% &
          interfaceEquationsToSolverMatrixMaps(interfaceMatrixIdx)%ptr,err,error,*999)
      ENDDO !interfaceMatrixIdx
      DEALLOCATE(interfaceToSolverMatrixMapsSm%interfaceEquationsToSolverMatrixMaps)
    ENDIF
    IF(ALLOCATED(interfaceToSolverMatrixMapsSm%interfaceColToSolverColsMap)) THEN
      DO columnIdx=1,SIZE(interfaceToSolverMatrixMapsSm%interfaceColToSolverColsMap,1)
        CALL SolverMapping_EquationsColToSolverColsMapFinalise(interfaceToSolverMatrixMapsSm% &
          & interfaceColToSolverColsMap(columnIdx),err,error,*999)
      ENDDO !columnIdx
      DEALLOCATE(interfaceToSolverMatrixMapsSm%interfaceColToSolverColsMap)
    ENDIF
   
    CALL Exits("SolverMapping_InterfaceToSolverMatrixMapsSmFinalise")
    RETURN
999 CALL Errors("SolverMapping_InterfaceToSolverMatrixMapsSmFinalise",err,error)    
    CALL Exits("SolverMapping_InterfaceToSolverMatrixMapsSmFinalise")
    RETURN 1
   
  END SUBROUTINE SolverMapping_InterfaceToSolverMatrixMapsSmFinalise

  !
  !================================================================================================================================
  !

  !>Initialises an interface to solver matrix maps sm.
  SUBROUTINE SolverMapping_InterfaceToSolverMatrixMapsSmInitialise(interfaceToSolverMatrixMapsSm,err,error,*)

    !Argument variables
    TYPE(InterfaceToSolverMatrixMapsSmType), INTENT(OUT) :: interfaceToSolverMatrixMapsSm !<The interface to solver matrix maps sm to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
                 
    CALL Enters("SolverMapping_InterfaceToSolverMatrixMapsSmInitialise",err,error,*999)

    interfaceToSolverMatrixMapsSm%solverMatrixNumber=0
    interfaceToSolverMatrixMapsSm%lagrangeVariableType=0
    NULLIFY(interfaceToSolverMatrixMapsSm%lagrangeVariable)
    CALL SolverMapping_VariableToSolverColMapInitialise(interfaceToSolverMatrixMapsSm%lagrangeVariableToSolverColMap,err,error,*999)
    interfaceToSolverMatrixMapsSm%numberOfDependentVariables=0
    interfaceToSolverMatrixMapsSm%numberOfInterfaceMatrices=0
        
    CALL Exits("SolverMapping_InterfaceToSolverMatrixMapsSmInitialise")
    RETURN
999 CALL Errors("SolverMapping_InterfaceToSolverMatrixMapsSmInitialise",err,error)    
    CALL Exits("SolverMapping_InterfaceToSolverMatrixMapsSmInitialise")
    RETURN 1
   
  END SUBROUTINE SolverMapping_InterfaceToSolverMatrixMapsSmInitialise

  !
  !================================================================================================================================
  !

  !>Finalises a Jacobian column to solver columns map and deallocates all memory.
  SUBROUTINE SolverMapping_JacobianColToSolverColsMapFinalise(jacobianColToSolverColsMap,err,error,*)

    !Argument variables
    TYPE(JacobianColToSolverColsMapType), INTENT(INOUT) :: jacobianColToSolverColsMap !<The Jacobian col to solver cols map to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    CALL Enters("SolverMapping_JacobianColToSolverColsMapFinalise",err,error,*999)

    IF(ALLOCATED(jacobianColToSolverColsMap%solverCols)) DEALLOCATE(jacobianColToSolverColsMap%solverCols)
    IF(ALLOCATED(jacobianColToSolverColsMap%couplingCoefficients)) DEALLOCATE(jacobianColToSolverColsMap%couplingCoefficients)
        
    CALL Exits("SolverMapping_JacobianColToSolverColsMapFinalise")
    RETURN
999 CALL Errors("SolverMapping_JacobianColToSolverColsMapFinalise",err,error)    
    CALL Exits("SolverMapping_JacobianColToSolverColsMapFinalise")
    RETURN 1
   
  END SUBROUTINE SolverMapping_JacobianColToSolverColsMapFinalise

  !
  !================================================================================================================================
  !

  !>Initialises an Jacobian column to solver columns map
  SUBROUTINE SolverMapping_JacobianColToSolverColsMapInitialise(jacobianColToSolverColsMap,err,error,*)

    !Argument variables
    TYPE(JacobianColToSolverColsMapType), INTENT(OUT) :: jacobianColToSolverColsMap !<The Jacobian column to solver columns map to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    CALL Enters("SolverMapping_JacobianColToSolverColsMapInitialise",err,error,*999)

    jacobianColToSolverColsMap%numberOfSolverCols=0
    
    CALL Exits("SolverMapping_JacobianColToSolverColsMapInitialise")
    RETURN
999 CALL Errors("SolverMapping_JacobianColToSolverColsMapInitialise",err,error)    
    CALL Exits("SolverMapping_JacobianColToSolverColsMapInitialise")
    RETURN 1
   
  END SUBROUTINE SolverMapping_JacobianColToSolverColsMapInitialise

  !
  !================================================================================================================================
  !

  !>Finalises a equations set to solver matrix map and deallocates all memory.
  SUBROUTINE SolverMapping_JacobianToSolverMapFinalise(jacobianToSolverMap,err,error,*)

    !Argument variables
    TYPE(JacobianToSolverMapType), POINTER :: jacobianToSolverMap !<The jacobian to solver map to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: columnIdx
    
    CALL Enters("SolverMapping_JacobianToSolverMapFinalise",err,error,*999)

    IF(ASSOCIATED(jacobianToSolverMap)) THEN
      IF(ALLOCATED(jacobianToSolverMap%jacobianColToSolverColsMap)) THEN
        DO columnIdx=1,SIZE(jacobianToSolverMap%jacobianColToSolverColsMap,1)
          CALL SolverMapping_JacobianColToSolverColsMapFinalise(jacobianToSolverMap% &
            & jacobianColToSolverColsMap(columnIdx),err,error,*999)
        ENDDO !columnIdx
        DEALLOCATE(jacobianToSolverMap%jacobianColToSolverColsMap)
      ENDIF
    ENDIF
        
    CALL Exits("SolverMapping_JacobianToSolverMapFinalise")
    RETURN
999 CALL Errors("SolverMapping_JacobianToSolverMapFinalise",err,error)    
    CALL Exits("SolverMapping_JacobianToSolverMapFinalise")
    RETURN 1
   
  END SUBROUTINE SolverMapping_JacobianToSolverMapFinalise

  !
  !================================================================================================================================
  !

  !>Initialises a Jacobian to solver maps
  SUBROUTINE SolverMapping_JacobianToSolverMapInitialise(jacobianToSolverMatrixMap,err,error,*)

    !Argument variables
    TYPE(JacobianToSolverMapType), POINTER :: jacobianToSolverMatrixMap !<The Jacobian to solver maps to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    CALL Enters("SolverMapping_JacobianToSolverMapInitialise",err,error,*999)

    IF(ASSOCIATED(jacobianToSolverMatrixMap)) THEN
      jacobianToSolverMatrixMap%solverMatrixNumber=0
      NULLIFY(jacobianToSolverMatrixMap%jacobianMatrix)
      NULLIFY(jacobianToSolverMatrixMap%solverMatrix)
    ELSE
      CALL FlagError("Jacobian to solver matrix map is not associated.",err,error,*999)
    ENDIF
        
    CALL Exits("SolverMapping_JacobianToSolverMapInitialise")
    RETURN
999 CALL Errors("SolverMapping_JacobianToSolverMapInitialise",err,error)    
    CALL Exits("SolverMapping_JacobianToSolverMapInitialise")
    RETURN 1
   
  END SUBROUTINE SolverMapping_JacobianToSolverMapInitialise

  !
  !================================================================================================================================
  !

  !>Finalises the solver column to dynamic equations map and deallocates all memory.
  SUBROUTINE SolverMapping_SolverColToDynamicEquationsMapFinalise(solverColToDynamicEquationsMap,err,error,*)

    !Argument variables
    TYPE(SolverColToDynamicEquationsMapType), INTENT(INOUT) :: solverColToDynamicEquationsMap !<The solver column to dynamic equations map to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    CALL Enters("SolverMapping_SolverColToDynamicEquationsMapFinalise",err,error,*999)

    IF(ALLOCATED(solverColToDynamicEquationsMap%equationsMatrixNumberS)) &
      & DEALLOCATE(solverColToDynamicEquationsMap%equationsMatrixNumbers)
    IF(ALLOCATED(solverColToDynamicEquationsMap%equationsColNumbers)) &
      & DEALLOCATE(solverColToDynamicEquationsMap%equationsColNumbers)
    IF(ALLOCATED(solverColToDynamicEquationsMap%couplingCoefficients)) &
      & DEALLOCATE(solverColToDynamicEquationsMap%couplingCoefficients)
       
    CALL Exits("SolverMapping_SolverColToDynamicEquationsMapFinalise")
    RETURN
999 CALL Errors("SolverMapping_SolverColToDynamicEquationsMapFinalise",err,error)
    CALL Exits("SolverMapping_SolverColToDynamicEquationsMapFinalise")
    RETURN 1
  END SUBROUTINE SolverMapping_SolverColToDynamicEquationsMapFinalise

  !
  !================================================================================================================================
  !

  !>Initialises the solver column to dynamic equations mapping and deallocates all memory.
  SUBROUTINE SolverMapping_SolverColToDynamicEquationsMapInitialise(solverColToDynamicEquationsMap,err,error,*)

    !Argument variables
    TYPE(SolverColToDynamicEquationsMapType), INTENT(OUT) :: solverColToDynamicEquationsMap !<The solver column to dynamic equations map to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
                 
    CALL Enters("SolverMapping_SolverColToDynamicEquationsMapInitialise",err,error,*999)

    solverColToDynamicEquationsMap%numberOfDynamicEquationsMatrices=0
     
    CALL Exits("SolverMapping_SolverColToDynamicEquationsMapInitialise")
    RETURN
999 CALL Errors("SolverMapping_SolverColToDynamicEquationsMapInitialise",err,error)
    CALL Exits("SolverMapping_SolverColToDynamicEquationsMapInitialise")
    RETURN 1
  END SUBROUTINE SolverMapping_SolverColToDynamicEquationsMapInitialise

  !
  !================================================================================================================================
  !

  !>Finalises the solver column to static equations map and deallocates all memory.
  SUBROUTINE SolverMapping_SolverColToStaticEquationsMapFinalise(solverColToStaticEquationsMap,err,error,*)

    !Argument variables
    TYPE(SolverColToStaticEquationsMapType), INTENT(INOUT) :: solverColToStaticEquationsMap !<The solver column to static equations map to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    CALL Enters("SolverMapping_SolverColToStaticEquationsMapFinalise",err,error,*999)

    IF(ALLOCATED(solverColToStaticEquationsMap%equationsMatrixNumbers)) &
      & DEALLOCATE(solverColToStaticEquationsMap%equationsMatrixNumbers)
    IF(ALLOCATED(solverColToStaticEquationsMap%equationsColNumbers)) &
      & DEALLOCATE(solverColToStaticEquationsMap%equationsColNumbers)
    IF(ALLOCATED(solverColToStaticEquationsMap%couplingCoefficients)) &
      & DEALLOCATE(solverColToStaticEquationsMap%couplingCoefficients)
       
    CALL Exits("SolverMapping_SolverColToStaticEquationsMapFinalise")
    RETURN
999 CALL Errors("SolverMapping_SolverColToStaticEquationsMapFinalise",err,error)
    CALL Exits("SolverMapping_SolverColToStaticEquationsMapFinalise")
    RETURN 1
  END SUBROUTINE SolverMapping_SolverColToStaticEquationsMapFinalise

  !
  !================================================================================================================================
  !

  !>Initialises the solver column to static equations mapping and deallocates all memory.
  SUBROUTINE SolverMapping_SolverColToStaticEquationsMapInitialise(solverColToStaticEquationsMap,err,error,*)

    !Argument variables
    TYPE(SolverColToStaticEquationsMapType), INTENT(OUT) :: solverColToStaticEquationsMap !<The solver column to static equations map to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    CALL Enters("SolverMapping_SolverColToStaticEquationsMapInitialise",err,error,*999)

    solverColToStaticEquationsMap%numberOfLinearEquationsMatrices=0
    solverColToStaticEquationsMap%jacobianColNumber=0    
    solverColToStaticEquationsMap%jacobianCouplingCoefficient=0.0_DP
    
    CALL Exits("SolverMapping_SolverColToStaticEquationsMapInitialise")
    RETURN
999 CALL Errors("SolverMapping_SolverColToStaticEquationsMapInitialise",err,error)
    CALL Exits("SolverMapping_SolverColToStaticEquationsMapInitialise")
    RETURN 1
  END SUBROUTINE SolverMapping_SolverColToStaticEquationsMapInitialise

  !
  !================================================================================================================================
  !

  !>Finalises the solver column to equations set map and deallocates all memory.
  SUBROUTINE SolverMapping_SolverColToEquationsSetMapFinalise(solverColToEquationsSetMap,err,error,*)

    !Argument variables
    TYPE(SolverColToEquationsSetMapType), INTENT(INOUT) :: solverColToEquationsSetMap !<The solver column to equations set map to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: columnIdx
    
    CALL Enters("SolverMapping_SolverColToEquationsSetMapFinalise",err,error,*999)

    IF(ALLOCATED(solverColToEquationsSetMap%solverColToDynamicEquationsMaps)) THEN
      DO columnIdx=1,SIZE(solverColToEquationsSetMap%solverColToDynamicEquationsMaps,1)
        CALL SolverMapping_SolverColToDynamicEquationsMapFinalise(solverColToEquationsSetMap% &
          & solverColToDynamicEquationsMaps(columnIdx),err,error,*999)
      ENDDO !columnIdx
      DEALLOCATE(solverColToEquationsSetMap%solverColToDynamicEquationsMaps)
      solverColToEquationsSetMap%haveDynamic=.FALSE.
    ENDIF
    IF(ALLOCATED(solverColToEquationsSetMap%SolverColToStaticEquationsMaps)) THEN
      DO columnIdx=1,SIZE(solverColToEquationsSetMap%SolverColToStaticEquationsMaps,1)
        CALL SolverMapping_SolverColToStaticEquationsMapFinalise(solverColToEquationsSetMap% &
          & SolverColToStaticEquationsMaps(columnIdx),err,error,*999)
      ENDDO !columnIdx
      DEALLOCATE(solverColToEquationsSetMap%SolverColToStaticEquationsMaps)
      solverColToEquationsSetMap%haveStatic=.FALSE.
    ENDIF
       
    CALL Exits("SolverMapping_SolverColToEquationsSetMapFinalise")
    RETURN
999 CALL Errors("SolverMapping_SolverColToEquationsSetMapFinalise",err,error)
    CALL Exits("SolverMapping_SolverColToEquationsSetMapFinalise")
    RETURN 1
  END SUBROUTINE SolverMapping_SolverColToEquationsSetMapFinalise

  !
  !================================================================================================================================
  !

  !>Initialises the solver column to equations set mapping and deallocates all memory.
  SUBROUTINE SolverMapping_SolverColToEquationsSetMapInitialise(solverColToEquationsSetMap,err,error,*)

    !Argument variables
    TYPE(SolverColToEquationsSetMapType), INTENT(OUT) :: solverColToEquationsSetMap !<The solver column to equations set map to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    CALL Enters("SolverMapping_SolverColToEquationsSetMapInitialise",err,error,*999)

    NULLIFY(solverColToEquationsSetMap%equations)
    solverColToEquationsSetMap%haveDynamic=.FALSE.
    solverColToEquationsSetMap%haveStatic=.FALSE.
    
    CALL Exits("SolverMapping_SolverColToEquationsSetMapInitialise")
    RETURN
999 CALL Errors("SolverMapping_SolverColToEquationsSetMapInitialise",err,error)
    CALL Exits("SolverMapping_SolverColToEquationsSetMapInitialise")
    RETURN 1
    
  END SUBROUTINE SolverMapping_SolverColToEquationsSetMapInitialise

  !
  !================================================================================================================================
  !

  !>Finalises the solver column to equations map and deallocates all memory.
  SUBROUTINE SolverMapping_SolverColToEquationsMapsFinalise(solverColToEquationsMaps,err,error,*)

    !Argument variables
    TYPE(SolverColToEquationsMapsType), INTENT(INOUT) :: solverColToEquationsMaps !<The solver column to equations sets map to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: columnIdx,equationsSetIdx
    
    CALL Enters("SolverMapping_SolverColToEquationsMapsFinalise",err,error,*999)

    IF(ALLOCATED(solverColToEquationsMaps%solverColToEquationsSetMaps)) THEN
      DO equationsSetIdx=1,SIZE(solverColToEquationsMaps%solverColToEquationsSetMaps,1)
        CALL SolverMapping_SolverColToEquationsSetMapFinalise(solverColToEquationsMaps% &
          & solverColToEquationsSetMaps(equationsSetIdx),err,error,*999)
      ENDDO !equationsSetIdx
      DEALLOCATE(solverColToEquationsMaps%solverColToEquationsSetMaps)
    ENDIF
    IF(ALLOCATED(solverColToEquationsMaps%solverDofToVariableMaps)) THEN
      DO columnIdx=1,SIZE(solverColToEquationsMaps%solverDofToVariableMaps,1)
        CALL SolverMapping_SolverDofToVariableMapFinalise(solverColToEquationsMaps%solverDofToVariableMaps(columnIdx), &
          & err,error,*999)
      ENDDO !columnIdx
      DEALLOCATE(solverColToEquationsMaps%solverDofToVariableMaps)
    ENDIF
    CALL DOMAIN_MAPPINGS_MAPPING_FINALISE(solverColToEquationsMaps%columnDofsMapping,err,error,*999)
    
    CALL Exits("SolverMapping_SolverColToEquationsMapsFinalise")
    RETURN
999 CALL Errors("SolverMapping_SolverColToEquationsMapsFinalise",err,error)
    CALL Exits("SolverMapping_SolverColToEquationsMapsFinalise")
    RETURN 1
  END SUBROUTINE SolverMapping_SolverColToEquationsMapsFinalise

  !
  !================================================================================================================================
  !

  !>Initialises the solver column to equations mapping and deallocates all memory.
  SUBROUTINE SolverMapping_SolverColToEquationsMapsInitialise(solverColToEquationsMaps,err,error,*)

    !Argument variables
    TYPE(SolverColToEquationsMapsType), INTENT(OUT) :: solverColToEquationsMaps !<The solver column to equations map to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
                 
    CALL Enters("SolverMapping_SolverColToEquationsMapsInitialise",err,error,*999)

    solverColToEquationsMaps%solverMatrixNumber=0
    NULLIFY(solverColToEquationsMaps%solverMatrix)
    NULLIFY(solverColToEquationsMaps%solverMapping)
    solverColToEquationsMaps%numberOfColumns=0
    solverColToEquationsMaps%numberOfDofs=0
    solverColToEquationsMaps%totalNumberOfDofs=0
    solverColToEquationsMaps%numberOfGlobalDofs=0
    NULLIFY(solverColToEquationsMaps%columnDofsMapping)
    
    CALL Exits("SolverMapping_SolverColToEquationsMapsInitialise")
    RETURN
999 CALL Errors("SolverMapping_SolverColToEquationsMapsiIntialise",err,error)
    CALL Exits("SolverMapping_SolverColToEquationsMapsiIntialise")
    RETURN 1
  END SUBROUTINE SolverMapping_SolverColToEquationsMapsInitialise

  !
  !================================================================================================================================
  !

  !>Finalises the solver column to interface map and deallocates all memory.
  SUBROUTINE SolverMapping_SolverColToInterfaceMapFinalise(solverColToInterfaceMap,err,error,*)

    !Argument variables
    TYPE(SolverColToInterfaceMapType), INTENT(INOUT) :: solverColToInterfaceMap !<The solver column to interface map to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: columnIdx
    
    CALL Enters("SolverMapping_SolverColToInterfaceMapFinalise",err,error,*999)

    NULLIFY(solverColToInterfaceMap%interfaceEquations)
    IF(ALLOCATED(solverColToInterfaceMap%solverColToInterfaceEquationsMaps)) THEN
      DO columnIdx=1,SIZE(solverColToInterfaceMap%solverColToInterfaceEquationsMaps,1)
        CALL SolverMapping_SolverColToInterfaceEquationsMapFinalise(solverColToInterfaceMap% &
          & solverColToInterfaceEquationsMaps(columnIdx),err,error,*999)
      ENDDO !columnIdx
      DEALLOCATE(solverColToInterfaceMap%solverColToInterfaceEquationsMaps)
    ENDIF
     
    CALL Exits("SolverMapping_SolverColToInterfaceMapFinalise")
    RETURN
999 CALL Errors("SolverMapping_SolverColToInterfaceMapFinalise",err,error)
    CALL Exits("SolverMapping_SolverColToInterfaceMapFinalise")
    RETURN 1
  END SUBROUTINE SolverMapping_SolverColToInterfaceMapFinalise

  !
  !================================================================================================================================
  !

  !>Initialises the solver column to interface mapping.
  SUBROUTINE SolverMapping_SolverColToInterfaceMapInitialise(solverColToInterfaceMap,err,error,*)

    !Argument variables
    TYPE(SolverColToInterfaceMapType), INTENT(OUT) :: solverColToInterfaceMap !<The solver column to interface map to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    CALL Enters("SolverMapping_SolverColToInterfaceMapInitialise",err,error,*999)

    NULLIFY(solverColToInterfaceMap%interfaceEquations)
    
    CALL Exits("SolverMapping_SolverColToInterfaceMapInitialise")
    RETURN
999 CALL Errors("SolverMapping_SolverColToInterfaceMapInitialise",err,error)
    CALL Exits("SolverMapping_SolverColToInterfaceMapInitialise")
    RETURN 1
    
  END SUBROUTINE SolverMapping_SolverColToInterfaceMapInitialise

  !
  !================================================================================================================================
  !

  !>Finalises the solver column to interface equations map and deallocates all memory.
  SUBROUTINE SolverMapping_SolverColToInterfaceEquationsMapFinalise(solverColToInterfaceEquationsMap,err,error,*)

    !Argument variables
    TYPE(SolverColToInterfaceEquationsMapType), INTENT(INOUT) :: solverColToInterfaceEquationsMap !<The solver column to interface equatiosn map to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    CALL Enters("SolverMapping_SolverColToInterfaceEquationsMapFinalise",err,error,*999)

    solverColToInterfaceEquationsMap%numberOfInterfaceMatrices=0
    IF(ALLOCATED(solverColToInterfaceEquationsMap%interfaceMatrixNumbers)) &
      & DEALLOCATE(solverColToInterfaceEquationsMap%interfaceMatrixNumbers)
    IF(ALLOCATED(solverColToInterfaceEquationsMap%interfaceColNumbers))  &
      & DEALLOCATE(solverColToInterfaceEquationsMap%interfaceColNumbers)
    IF(ALLOCATED(solverColToInterfaceEquationsMap%couplingCoefficients))  &
      & DEALLOCATE(solverColToInterfaceEquationsMap%couplingCoefficients)
    
    CALL Exits("SolverMapping_SolverColToInterfaceEquationsMapFinalise")
    RETURN
999 CALL Errors("SolverMapping_SolverColToInterfaceEquationsMapFinalise",err,error)
    CALL Exits("SolverMapping_SolverColToInterfaceEquationsMapFinalise")
    RETURN 1
  END SUBROUTINE SolverMapping_SolverColToInterfaceEquationsMapFinalise

  !
  !================================================================================================================================
  !

  !>Initialises the solver column to interface equations mapping.
  SUBROUTINE SolverMapping_SolverColToInterfaceEquationsMapInitialise(solverColToInterfaceEquationsMap,err,error,*)

    !Argument variables
    TYPE(SolverColToInterfaceEquationsMapType), INTENT(OUT) :: solverColToInterfaceEquationsMap !<The solver column to interface equations map to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    CALL Enters("SolverMapping_SolverColToInterfaceEquationsMapInitialise",err,error,*999)

    solverColToInterfaceEquationsMap%numberOfInterfaceMatrices=0
    
    CALL Exits("SolverMapping_SolverColToInterfaceEquationsMapInitialise")
    RETURN
999 CALL Errors("SolverMapping_SolverColToInterfaceEquationsMapInitialise",err,error)
    CALL Exits("SolverMapping_SolverColToInterfaceEquationsMapInitialise")
    RETURN 1
    
  END SUBROUTINE SolverMapping_SolverColToInterfaceEquationsMapInitialise

  !
  !================================================================================================================================
  !

  !>Finalises the solver dof to variable mapping and deallocates all memory.
  SUBROUTINE SolverMapping_SolverDofToVariableMapFinalise(solverDofToVariableMap,err,error,*)

    !Argument variables
    TYPE(SolverDofToVariableMapType) :: solverDofToVariableMap !<The solver dof to variable map to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    CALL Enters("SolverMapping_SolverDofToVariableMapFinalise",err,error,*999)

    IF(ALLOCATED(solverDofToVariableMap%equationsTypes)) DEALLOCATE(solverDofToVariableMap%equationsTypes)
    IF(ALLOCATED(solverDofToVariableMap%equationsIndices)) DEALLOCATE(solverDofToVariableMap%equationsIndices)
    IF(ALLOCATED(solverDofToVariableMap%VARIABLE)) DEALLOCATE(solverDofToVariableMap%VARIABLE)
    IF(ALLOCATED(solverDofToVariableMap%variableDof)) DEALLOCATE(solverDofToVariableMap%variableDof)
    IF(ALLOCATED(solverDofToVariableMap%variableCoefficient)) DEALLOCATE(solverDofToVariableMap%variableCoefficient)
    IF(ALLOCATED(solverDofToVariableMap%additiveConstant)) DEALLOCATE(solverDofToVariableMap%additiveConstant)
    
    CALL Exits("SolverMapping_SolverDofToVariableMapFinalise")
    RETURN
999 CALL Errors("SolverMapping_SolverDofToVariableMapFinalise",err,error)
    CALL Exits("SolverMapping_SolverDofToVariableMapFinalise")
    RETURN 1
  END SUBROUTINE SolverMapping_SolverDofToVariableMapFinalise

  !
  !================================================================================================================================
  !

  !>Initialises the solver dof to variable mapping.
  SUBROUTINE SolverMapping_SolverDofToVariableMapInitialise(solverDofToVariableMap,err,error,*)

    !Argument variables
    TYPE(SolverDofToVariableMapType) :: solverDofToVariableMap !<The solver dof to variable map to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    CALL Enters("SolverMapping_SolverDofToVariableMapInitialise",err,error,*999)

    solverDofToVariableMap%numberOfEquations=0
    
    CALL Exits("SolverMapping_SolverDofToVariableMapInitialise")
    RETURN
999 CALL Errors("SolverMapping_SolverDofToVariableMapInitialise",err,error)
    CALL Exits("SolverMapping_SolverDofToVariableMapInitialise")
    RETURN 1
  END SUBROUTINE SolverMapping_SolverDofToVariableMapInitialise

  !
  !================================================================================================================================
  !

  !>Sets/changes the number of solver matrices for the solver mapping
  SUBROUTINE SolverMapping_SolverMatricesNumberSet(solverMapping,numberOfSolverMatrices,err,error,*)

    !Argument variables
    TYPE(SolverMappingType), POINTER :: solverMapping !<A pointer to the solver mapping
    INTEGER(INTG), INTENT(IN) :: numberOfSolverMatrices !<The number of solver matrices for the solver.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: equationsSetIdx,matrixIdx,maximumNumberOfEquationsMatrices
    INTEGER(INTG), ALLOCATABLE :: newMatrixVariableTypes(:,:,:)
    TYPE(EQUATIONS_TYPE), POINTER :: equations
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: equationsMapping
    TYPE(EQUATIONS_MAPPING_LINEAR_TYPE), POINTER :: linearMapping
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet
    TYPE(VARYING_STRING) :: localError

    CALL Enters("SolverMapping_SolverMatricesNumberSet",err,error,*999)

    IF(ASSOCIATED(solverMapping)) THEN
      IF(solverMapping%solverMappingFinished) THEN
        CALL FlagError("Solver mappings has been finished.",err,error,*999)
      ELSE
        IF(ASSOCIATED(solverMapping%createValuesCache)) THEN
          maximumNumberOfEquationsMatrices=1
          DO equationsSetIdx=1,solverMapping%numberOfEquationsSets
            equationsSet=>solverMapping%equationsSets(equationsSetIdx)%ptr
            IF(ASSOCIATED(equationsSet)) THEN
              equations=>equationsSet%equations
              IF(ASSOCIATED(equations)) THEN
                equationsMapping=>equations%EQUATIONS_MAPPING
                IF(ASSOCIATED(equationsMapping)) THEN
                  linearMapping=>equationsMapping%LINEAR_MAPPING
                  IF(ASSOCIATED(linearMapping)) THEN
                    IF(ASSOCIATED(linearMapping)) THEN
                      IF(linearMapping%NUMBER_OF_LINEAR_EQUATIONS_MATRICES>maximumNumberOfEquationsMatrices) &
                        & maximumNumberOfEquationsMatrices=linearMapping%NUMBER_OF_LINEAR_EQUATIONS_MATRICES
                    ENDIF
                  ELSE
                    CALL FlagError("Equations mapping linear mapping is not associated.",err,error,*999)
                  ENDIF
                ELSE
                  CALL FlagError("Equations equations mapping is not associated.",err,error,*999)
                ENDIF
              ELSE
                CALL FlagError("Equations set equations is not associated.",err,error,*999)
              ENDIF
            ELSE
              localError="Equations set is not associated for equations set index "// &
                & TRIM(NumberToVString(equationsSetIdx,"*",err,error))//"."
              CALL FlagError(localError,err,error,*999)
            ENDIF
          ENDDO !equationsSetIdx
          !Check number of matrices to set is valid
          IF(numberOfSolverMatrices>0.AND.numberOfSolverMatrices<=maximumNumberOfEquationsMatrices) THEN
            !If we need to reallocate and reset all the create values cache arrays and change the number of matrices
            IF(numberOfSolverMatrices/=solverMapping%numberOfSolverMatrices) THEN
              ALLOCATE(newMatrixVariableTypes(0:FIELD_NUMBER_OF_VARIABLE_TYPES,solverMapping%numberOfEquationsSets, &
                & numberOfSolverMatrices),STAT=err)
              IF(err/=0) CALL FlagError("Could not allocate old matrix variable types",err,error,*999)              
              IF(numberOfSolverMatrices>solverMapping%numberOfSolverMatrices) THEN
                newMatrixVariableTypes(:,:,1:solverMapping%numberOfSolverMatrices)=solverMapping%createValuesCache% &
                  & matrixVariableTypes(:,:,1:solverMapping%numberOfSolverMatrices)
                DO matrixIdx=solverMapping%numberOfSolverMatrices+1,numberOfSolverMatrices
                  newMatrixVariableTypes(0,:,matrixIdx)=1
                  newMatrixVariableTypes(1,:,matrixIdx)=FIELD_U_VARIABLE_TYPE
                  newMatrixVariableTypes(2:FIELD_NUMBER_OF_VARIABLE_TYPES,:,matrixIdx)=0
                ENDDO !matrixIdx
              ELSE
                newMatrixVariableTypes(:,:,1:numberOfSolverMatrices)=solverMapping%createValuesCache% &
                  & matrixVariableTypes(:,:,1:numberOfSolverMatrices)
              ENDIF
              CALL MOVE_ALLOC(newMatrixVariableTypes,solverMapping%createValuesCache%matrixVariableTypes)
              solverMapping%numberOfSolverMatrices=numberOfSolverMatrices
            ENDIF
          ELSE
            localError="The specified number of solver matrices of "// &
              & TRIM(NumberToVString(numberOfSolverMatrices,"*",err,error))// &
              & " is invalid. The number must be >= 1 and <= "// &
              & TRIM(NumberToVString(maximumNumberOfEquationsMatrices,"*",err,error))//"."
            CALL FlagError(localError,err,error,*999)
          ENDIF
        ELSE
          CALL FlagError("Solver mapping create values cache is not associated",err,error,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Solver mapping is not associated",err,error,*999)
    ENDIF
    
    CALL Exits("SolverMapping_SolverMatricesNumberSet")
    RETURN
999 IF(ALLOCATED(newMatrixVariableTypes)) DEALLOCATE(newMatrixVariableTypes)
    CALL Errors("SolverMapping_SolverMatricesNumberSet",err,error)
    CALL Exits("SolverMapping_SolverMatricesNumberSet")
    RETURN 1
  END SUBROUTINE SolverMapping_SolverMatricesNumberSet
  
  !
  !================================================================================================================================
  !

  !>Finalises a solver row to equations map and deallocates all memory.
  SUBROUTINE SolverMapping_SolverRowToEquationsMapsFinalise(solverRowToEquationsMaps,err,error,*)

    !Argument variables
    TYPE(SolverRowToEquationsMapsType) :: solverRowToEquationsMaps !<The solver row to equations map to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    CALL Enters("SolverMapping_SolverRowToEquationsMapsFinalise",err,error,*999)

    solverRowToEquationsMaps%numberOfEquationsSets=0
    solverRowToEquationsMaps%interfaceConditionIndex=0
    IF(ALLOCATED(solverRowToEquationsMaps%equationsIndex)) DEALLOCATE(solverRowToEquationsMaps%equationsIndex)
    IF(ALLOCATED(solverRowToEquationsMaps%rowcolNumber)) DEALLOCATE(solverRowToEquationsMaps%rowcolNumber)
    IF(ALLOCATED(solverRowToEquationsMaps%couplingCoefficients)) DEALLOCATE(solverRowToEquationsMaps%couplingCoefficients)
    
    CALL Exits("SolverMapping_SolverRowToEquationsMapsFinalise")
    RETURN
999 CALL Errors("SolverMapping_SolverRowToEquationsMapsFinalise",err,error)    
    CALL Exits("SolverMapping_SolverRowToEquationsMapsFinalise")
    RETURN 1
    
  END SUBROUTINE SolverMapping_SolverRowToEquationsMapsFinalise

  !
  !================================================================================================================================
  !

  !>Initialises a solver row to equations map.
  SUBROUTINE SolverMapping_SolverRowToEquationsMapsInitialise(solverRowToEquationsMaps,err,error,*)

    !Argument variables
    TYPE(SolverRowToEquationsMapsType), INTENT(OUT) :: solverRowToEquationsMaps !<The solver row to equations map to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    CALL Enters("SolverMapping_SolverRowToEquationsMapsInitialise",err,error,*999)

    solverRowToEquationsMaps%numberOfEquationsSets=0
    solverRowToEquationsMaps%interfaceConditionIndex=0
        
    CALL Exits("SolverMapping_SolverRowToEquationsMapsInitialise")
    RETURN
999 CALL Errors("SolverMapping_SolverRowToEquationsMapsInitialise",err,error)    
    CALL Exits("SolverMapping_SolverRowToEquationsMapsInitialise")
    RETURN 1
    
  END SUBROUTINE SolverMapping_SolverRowToEquationsMapsInitialise

  !
  !================================================================================================================================
  !

  !>Adds a field variable to a list of solver variables
  SUBROUTINE SolverMapping_SolverVariableAdd(solverVariables,variable,equationsType,equationsIdx,err,error,*)

    !Argument variables
    TYPE(SolverMappingVariablesType), INTENT(INOUT) :: solverVariables !<The solver variables to add the variable to
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: variable !<The field varaible to add.
    INTEGER(INTG), INTENT(IN) :: equationsType !<The type of equations to add \see SolverMapping_EquationsTypes,SolverMapping
    INTEGER(INTG), INTENT(IN) :: equationsIdx !<The index of the equations to add
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: solverVariableIdx
    TYPE(SolverMappingVariablePtrType), ALLOCATABLE :: newSolverVariables
 
    CALL Enters("SolverMapping_SolverVariableAdd",err,error,*997)

    IF(ASSOCIATED(variable)) THEN
      ALLOCATE(newSolverVariables(solverVariables%numberOfVariables+1),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocated new solver variables.",err,error,*999)
      DO solverVariableIdx=1,solverVariables%numberOfVariables
        newSolverVariables(solverVariableIdx)%ptr=>solverVariables%variables(solverVariableIdx)%ptr
      ENDDO !solverVariableIdx
      NULLIFY(newSolverVariables(solverVariables%numberOfVariables+1)%ptr)
      CALL SolverMapping_SolverVariableInitialise(newSolverVariables(solverVariables%numberOfVariables+1)%ptr,err,error,*999)
      newSolverVariables(solverVariables%numberOfVariables+1)%ptr%variable=>variable
      newSolverVariables(solverVariables%numberOfVariables+1)%ptr%variableType=>variable%TYPE
      CALL SolverMapping_SolverVariableEquationsAdd(newSolverVariables(solverVariables%numberOfVariables+1)%ptr, &
        & equationsType,equationsIdx,err,error,*999)
      CALL MOVE_ALLOC(newSolverVariables,solverVariables%variables)
      solverVariables%numberOfVariables=solverVariables%numberOfVariables+1
    ELSE
      CALL FlagError("Variable is not associated.",err,error,*999)
    ENDIF
    
    CALL Exits("SolverMapping_SolverVariableAdd")
    RETURN
999 IF(ALLOCATED(newSolverVariables)) THEN
      CALL SolverMapping_SolverVariableFinalise(newSolverVariables(solverVariables%numberOfVariables+1)%ptr,err,error,*998)
998   DEALLOCATE(newSolverVariables)
    ENDIF
997 CALL Errors("SolverMapping_SolverVariableAdd",err,error)
    CALL Exits("SolverMapping_SolverVariableAdd")
    RETURN 1
  END SUBROUTINE SolverMapping_SolverVariableAdd

  !
  !================================================================================================================================
  !

  !>Adds an equations to the list of equations involving a particular solver variable.
  SUBROUTINE SolverMapping_SolverVariableEquationsAdd(solverVariable,equationsType,equationsIdx,err,error,*)

    !Argument variables
    TYPE(SolverMappingVariableType), INTENT(INOUT) :: solverVariable !<The solver variables to add the equations to
    INTEGER(INTG), INTENT(IN) :: equationsType !<The type of equations to add \see SolverMapping_EquationsTypes,SolverMapping
    INTEGER(INTG), INTENT(IN) :: equationsIdx !<The index of the equations to add
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG), ALLOCATABLE :: newEquationIndices,newEquationTypes

    CALL Enters("SolverMapping_SolverVariableEquationsAdd",err,error,*999)

    ALLOCATE(newEquationTypes(solverVariable%numberOfEquations+1),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate new equation types.",err,error,*999)
    ALLOCATE(newEquationIndices(solverVariable%numberOfEquations+1),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate new equation indices.",err,error,*999)
    newEquationTypes(1:solverVariable%numberOfEquations)=solverVariable%equationTypes(1:solverVariable%numberOfEquations)
    newEquationIndices(1:solverVariable%numberOfEquations)=solverVariable%equationIndices(1:solverVariable%numberOfEquations)
    newEquationTypes(solverVariable%numberOfEquations+1)=equationsType
    newEquationIndices(solverVariable%numberOfEquations+1)=equationsIdx
    CALL MOVE_ALLOC(newEquationTypes,solverVariable%equationTypes)
    CALL MOVE_ALLOC(newEquationIndices,solverVariable%equationIndices)
    solverVariable%numberOfEquations=solverVariable%numberOfEquations+1
    
    CALL Exits("SolverMapping_SolverVariableEquationsAdd")
    RETURN
999 IF(ALLOCATED(newEquationTypes)) DEALLOCATE(newEquationTypes)
    IF(ALLOCATED(newEquationIndices)) DEALLOCATE(newEquationIndices)    
    CALL Errors("SolverMapping_SolverVariableEquationsAdd",err,error)
    CALL Exits("SolverMapping_SolverVariableEquationsAdd")
    RETURN 1
  END SUBROUTINE SolverMapping_SolverVariableEquationsAdd

  !
  !================================================================================================================================
  !

  !>Finalises the solver mapping variable and deallocates all memory.
  SUBROUTINE SolverMapping_VariableFinalise(solverVariable,err,error,*)

    !Argument variables
    TYPE(SolverMappingVariableType), POINTER :: solverVariable !<The solver mapping variable to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    CALL Enters("SolverMapping_VariableFinalise",err,error,*999)

    IF(ASSOCIATED(solverVariable)) THEN
      NULLIFY(solverMappingVariable%variable)
      solverMappingVariable%variableType=0
      solverMappingVariable%numberOfEquations=0
      IF(ALLOCATED(solverMappingVariable%equationTypes)) DEALLOCATE(solverMappingVariable%equationTypes)
      IF(ALLOCATED(solverMappingVariable%equationIndices)) DEALLOCATE(solverMappingVariable%equationIndices)
      DEALLOCATE(solverVariable)
    ENDIF
    
    CALL Exits("SolverMapping_VariableFinalise")
    RETURN
999 CALL Errors("SolverMapping_VariableFinalise",err,error)
    CALL Exits("SolverMapping_VariableFinalise")
    RETURN 1
  END SUBROUTINE SolverMapping_VariableFinalise

  !
  !================================================================================================================================
  !

  !>Trys to find a field variable in a list of solver variables. If it finds the variable the it adds the equations to the
  !>list of equations involving that solver variable. If it doesn't find the variable a new solver variable is added to the
  !>list of solver variables.
  SUBROUTINE SolverMapping_SolverVariableFindAdd(solverVariablesList,variable,equationsType,equationsIdx,err,error,*)

    !Argument variables
    TYPE(SolverMappingVariablesType), INTENT(INOUT) :: solverVariablesList !<The list of solver variables to find or add to
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: variable !<The field variable to find or add
    INTEGER(INTG), INTENT(IN) :: equationsType !<The type of equations the field variable is from \see SolverMapping_EquationsTypes,SolverMapping
    INTEGER(INTG), INTENT(IN) :: equationsIdx !<The index of the equations the field variable is from
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: solverVariableIdx
    LOGICAL :: found
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: solverVariable
    TYPE(VARYING_STRING) :: localError

    CALL Enters("SolverMapping_SolverVariableFindAdd",err,error,*999)

    IF(ASSOCIATED(variable)) THEN
      found=.FALSE.
      DO solverVariableIdx=1,solverVariablesList%numberOfVariables
        solverVariable=>solverVariablesList%variables(solverVariableIdx)%variable
        IF(ASSOCIATED(solverVariable)) THEN
          IF(ASSOCIATED(solverVariable,variable)) THEN
            !Found the variable, add the equations set to the list of equations involving the solver variable
            found=.TRUE.
            CALL SolverMapping_SolverVariableEquationsAdd(solverVariablesList%variables(solverVariableIdx), &
              & equationsType,equationsIdx,err,error,*999)
            EXIT
          ENDIF
        ELSE
          localError="Solver variable is not associated for solver variable index "// &
            & TRIM(NumberToVstring(variableIdx,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ENDDO !variableIdx
      IF(.NOT.found) THEN
        !New solver variable, add it to the list of row variables.
        CALL SolverMapping_SolverVariableAdd(solveVariablesList,variable,equationsType,equationsIdx,err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Variable is not associated.",err,error,*999)
    ENDIF
    
    CALL Exits("SolverMapping_SolverVariableFindAdd")
    RETURN
999 CALL Errors("SolverMapping_SolverVariableFindAdd",err,error)
    CALL Exits("SolverMapping_SolverVariableFindAdd")
    RETURN 1
  END SUBROUTINE SolverMapping_SolverVariableFindAdd

  !
  !================================================================================================================================
  !

  !>Initialises the solver mapping and deallocates all memory.
  SUBROUTINE SolverMapping_VariableInitialise(solverVariable,err,error,*)

    !Argument variables
    TYPE(SolverMappingVariableType), POINTER :: solverVariable !<The solver mapping variable to initialise. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError
 
    CALL Enters("SolverMapping_VariableInitialise",err,error,*998)

    IF(ASSOCIATED(solverVariable)) THEN
      CALL FlagError("Solver variable is already associated.",err,error,*998)
    ELSE
      ALLOCATE(solverVariable,STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate solver variable.",err,error,*999)
      NULLIFY(solverMappingVariable%variable)
      solverMappingVariable%variableType=0
      solverMappingVariable%numberOfEquations=0
    ENDIF
    
    CALL Exits("SolverMapping_VariableInitialise")
    RETURN
999 CALL SolverMapping_VariableFinalise(solverVariable,dummyErr,dummyError,*998)
998 CALL Errors("SolverMapping_VariableInitialise",err,error)
    CALL Exits("SolverMapping_VariableInitialise")
    RETURN 1
  END SUBROUTINE SolverMapping_VariableInitialise

  !
  !================================================================================================================================
  !

  !>Finalises the solver mapping variables and deallocates all memory.
  SUBROUTINE SolverMapping_VariablesFinalise(solverMappingVariables,err,error,*)

    !Argument variables
    TYPE(SolverMappingVariablesType), INTENT(INOUT) :: solverMappingVariables !<The solver mapping variables to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: variableIdx

    CALL Enters("SolverMapping_VariablesFinalise",err,error,*999)

    solverMappingVariables%numberOfVariables=0
    IF(ALLOCATED(solverMappingVariables%variables)) THEN
      DO variableIdx=1,SIZE(solverMappingVariables%variables,1)
        CALL SolverMapping_VariableFinalise(solverMappingVariables%variables(variableIdx)%ptr,err,error,*999)
      ENDDO !variableIdx
      DEALLOCATE(solverMappingVariables%variables)
    ENDIF
     
    CALL Exits("SolverMapping_VariablesFinalise")
    RETURN
999 CALL Errors("SolverMapping_VariablesFinalise",err,error)
    CALL Exits("SolverMapping_VariablesFinalise")
    RETURN 1
  END SUBROUTINE SolverMapping_VariablesFinalise

  !
  !================================================================================================================================
  !

  !>Initialises the solver mapping variables and deallocates all memory.
  SUBROUTINE SolverMapping_VariablesInitialise(solverMappingVariables,err,error,*)

    !Argument variables          !Allocate the column variable lists
    TYPE(SolverMappingVariablesType), INTENT(OUT) :: solverMappingVariables !<The solver mapping variables to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    CALL Enters("SolverMapping_VariablesInitialise",err,error,*999)

    solverMappingVariables%numberOfVariables=0
    
    CALL Exits("SolverMapping_VariablesInitialise")
    RETURN
999 CALL Errors("SolverMapping_VariablesInitialise",err,error)
    CALL Exits("SolverMapping_VariablesInitialise")
    RETURN 1
  END SUBROUTINE SolverMapping_VariablesInitialise

  !
  !================================================================================================================================
  !

  !>Finalises a variable to solver column map and deallocates all memory.
  SUBROUTINE SolverMapping_VariableToSolverColMapFinalise(variableToSolverColMap,err,error,*)

    !Argument variables
    TYPE(VariableToSolverColMapType), INTENT(INOUT) :: variableToSolverColMap !<The variable to solver column map to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL Enters("SolverMapping_VariableToSolverColMapFinalise",err,error,*999)

    IF(ALLOCATED(variableToSolverColMap%columnNumbers)) DEALLOCATE(variableToSolverColMap%columnNumbers)
    IF(ALLOCATED(variableToSolverColMap%couplingCoefficients)) DEALLOCATE(variableToSolverColMap%couplingCoefficients)
    IF(ALLOCATED(variableToSolverColMap%additiveConstants)) DEALLOCATE(variableToSolverColMap%additiveConstants)
        
    CALL Exits("SolverMapping_VariableToSolverColMapFinalise")
    RETURN
999 CALL Errors("SolverMapping_VariableToSolverColMapFinalise",err,error)    
    CALL Exits("SolverMapping_VariableToSolverColMapFinalise")
    RETURN 1
   
  END SUBROUTINE SolverMapping_VariableToSolverColMapFinalise

  !
  !================================================================================================================================
  !

  !>Initialises a variable to solver column map.
  SUBROUTINE SolverMapping_VariableToSolverColMapInitialise(variableToSolverColMap,err,error,*)

    !Argument variables
    TYPE(VariableToSolverColMapType), INTENT(IN) :: variableToSolverColMap  !<The variable to solver column map to initalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    CALL Enters("SolverMapping_VariableToSolverColMapInitialise",err,error,*999)

    !Nothing to do here, all members are allocatable

    CALL Exits("SolverMapping_VariableToSolverColMapInitialise")
    RETURN
999 CALL Errors("SolverMapping_VariableToSolverColMapInitialise",err,error)    
    CALL Exits("SolverMapping_VariableToSolverColMapInitialise")
    RETURN 1
   
  END SUBROUTINE SolverMapping_VariableToSolverColMapInitialise

  !
  !================================================================================================================================
  !

  !>Add a new DOF coupling to the list of global solver couplings
  SUBROUTINE SolverDofCouplings_AddCoupling(solverDofCouplings,coupling,couplingIndex,err,error,*)

    !Argument variables
    TYPE(SolverMappingDofCouplingsType), INTENT(INOUT) :: solverDofCouplings !<The solver row or column couplings
    TYPE(BoundaryConditionsCoupledDofsType), POINTER, INTENT(IN) :: coupling !<The new DOF boundary condition coupling to add to the list
    INTEGER(INTG), INTENT(OUT) :: couplingIndex !<On return, the index of the coupling in the solver dof couplings
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: newSize
    TYPE(BoundaryConditionsCoupledDofsPtrType), ALLOCATABLE :: newDofCouplings(:)

    CALL Enters("SolverDofCouplings_AddCoupling",err,error,*998)

    IF(solverDofCouplings%numberOfCouplings+1>solverDofCouplings%capacity) THEN
      !Resize or perform initial allocation if necessary
      IF(solverDofCouplings%capacity==0) THEN
        newSize=32
      ELSE
        newSize=solverDofCouplings%capacity*2
      END IF
      ALLOCATE(newDofCouplings(newSize),stat=err)
      IF(err/=0) CALL FlagError("Could not allocate new DOF couplings array.",err,error,*999)
      IF(solverDofCouplings%capacity>0) THEN
        newDofCouplings(1:solverDofCouplings%numberOfCouplings)= &
          & solverDofCouplings%dofCouplings(1:solverDofCouplings%numberOfCouplings)
      END IF
      solverDofCouplings%capacity=newSize
      CALL MOVE_ALLOC(newDofCouplings,solverDofCouplings%dofCouplings)
    END IF

    solverDofCouplings%dofCouplings(solverDofCouplings%numberOfCouplings+1)%ptr=>coupling
    solverDofCouplings%numberOfCouplings=solverDofCouplings%numberOfCouplings+1
    couplingIndex=solverDofCouplings%numberOfCouplings

    CALL Exits("SolverDofCouplings_AddCoupling")
    RETURN
999 IF(ALLOCATED(newDofCouplings)) DEALLOCATE(newDofCouplings)
998 CALL Errors("SolverDofCouplings_AddCoupling",err,error)
    CALL Exits("SolverDofCouplings_AddCoupling")
    RETURN 1

  END SUBROUTINE SolverDofCouplings_AddCoupling

  !
  !================================================================================================================================
  !

  !>Initialise the list of solver row or column couplings
  SUBROUTINE SolverDofCouplings_Finalise(solverDofCouplings,err,error,*)

    !Argument variables
    TYPE(SolverMappingDofCouplingsType), INTENT(INOUT) :: solverDofCouplings !<The solver row or column couplings to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string

    CALL Enters("SolverDofCouplings_Finalise",err,error,*999)

    IF(ALLOCATED(solverDofCouplings%dofCouplings)) THEN
      !Don't finalise individual DOF couplings as these are owned
      !by the boundary conditions.
      DEALLOCATE(solverDofCouplings%dofCouplings)
    END IF
    solverDofCouplings%numberOfCouplings=0
    solverDofCouplings%capacity=0

    CALL Exits("SolverDofCouplings_Finalise")
    RETURN
999 CALL Errors("SolverDofCouplings_Finalise",err,error)
    CALL Exits("SolverDofCouplings_Finalise")
    RETURN 1

  END SUBROUTINE SolverDofCouplings_Finalise

  !
  !================================================================================================================================
  !

  !>Initialise the list of solver row or column couplings
  SUBROUTINE SolverDofCouplings_Initialise(solverDofCouplings,err,error,*)

    !Argument variables
    TYPE(SolverMappingDofCouplingsType), INTENT(INOUT) :: solverDofCouplings !<The solver row or column couplings to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string

    CALL Enters("SolverDofCouplings_Initialise",err,error,*999)

    solverDofCouplings%numberOfCouplings=0
    solverDofCouplings%capacity=0

    CALL Exits("SolverDofCouplings_Initialise")
    RETURN
999 CALL Errors("SolverDofCouplings_Initialise",err,error)
    CALL Exits("SolverDofCouplings_Initialise")
    RETURN 1

  END SUBROUTINE SolverDofCouplings_Initialise

END MODULE SolverMappingRoutines
