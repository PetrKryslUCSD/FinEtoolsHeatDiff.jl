using SparseArrays
import FinEtools.AssemblyModule: AbstractSysmatAssembler, startassembly!, assemble!, makematrix!
using Sparspak
using Sparspak.SpkProblem: Problem, inaij!
using Sparspak.SpkSparseSolver: SparseSolver, inmatrix!, findorder!, symbolicfactor!
using Sparspak.SparseSolver: SparseSolver, solve!
using MKL

"""
    SysmatAssemblerSparspak{IT, MBT, IBT} <: AbstractSysmatAssembler

Type for assembling a sparse global matrix from elementwise matrices.

!!! note

    All fields of the datatype are private. The type is manipulated by the
    functions `startassembly!`, `assemble!`, and `makematrix!`.
"""
mutable struct SysmatAssemblerSparspak{IT, FT} <: AbstractSysmatAssembler
    problem::Problem{IT, FT}
    ndofs_row::IT
    ndofs_col::IT
    nomatrixresult::Bool
end

"""
    SysmatAssemblerSparspak(z= zero(FFlt), nomatrixresult = false)

Construct blank system matrix assembler.

The matrix entries are of type `T`. The assembler either produces a sparse
matrix (when `nomatrixresult = true`), or does not (when `nomatrixresult =
false`). When the assembler does not produce the sparse matrix when
`makematrix!` is called, it still can be constructed from the buffers stored in
the assembler.


# Example

This is how a sparse matrix is assembled from two rectangular dense matrices.
```
    a = SysmatAssemblerSparspak(0.0)
    startassembly!(a, 5, 5, 3, 7, 7)
    m = [0.24406   0.599773    0.833404  0.0420141
        0.786024  0.00206713  0.995379  0.780298
        0.845816  0.198459    0.355149  0.224996]
    assemble!(a, m, [1 7 5], [5 2 1 4])
    m = [0.146618  0.53471   0.614342    0.737833
         0.479719  0.41354   0.00760941  0.836455
         0.254868  0.476189  0.460794    0.00919633
         0.159064  0.261821  0.317078    0.77646
         0.643538  0.429817  0.59788     0.958909]
    assemble!(a, m, [2 3 1 7 5], [6 7 3 4])
    A = makematrix!(a)
```

When the `nomatrixresult` is set as true, no matrix is produced.
```
    a = SysmatAssemblerSparspak(0.0, true)
    startassembly!(a, 5, 5, 3, 7, 7)
    m = [0.24406   0.599773    0.833404  0.0420141
        0.786024  0.00206713  0.995379  0.780298
        0.845816  0.198459    0.355149  0.224996]
    assemble!(a, m, [1 7 5], [5 2 1 4])
    m = [0.146618  0.53471   0.614342    0.737833
         0.479719  0.41354   0.00760941  0.836455
         0.254868  0.476189  0.460794    0.00919633
         0.159064  0.261821  0.317078    0.77646
         0.643538  0.429817  0.59788     0.958909]
    assemble!(a, m, [2 3 1 7 5], [6 7 3 4])
    A = makematrix!(a)
```
Here `A` is a sparse zero matrix. To construct the correct matrix is still
possible, for instance like this:
```
    a.nomatrixresult = false
    A = makematrix!(a)
```
At this point all the buffers of the assembler have been cleared, and
`makematrix!(a) ` is no longer possible.

"""
function SysmatAssemblerSparspak(z= zero(FFlt), nomatrixresult = false)
    return SysmatAssemblerSparspak(Problem(0, 0, 2500, z), 0, 0, nomatrixresult)
end

"""
    startassembly!(
        self::SysmatAssemblerSparspak,
        elem_mat_nrows,
        elem_mat_ncols,
        elem_mat_nmatrices,
        ndofs_row,
        ndofs_col,
    )

Start the assembly of a global matrix.

The method makes buffers for matrix assembly. It must be called before
the first call to the method `assemble!`.
- `elem_mat_nrows`= number of rows in typical element matrix,
- `elem_mat_ncols`= number of columns in a typical element matrix,
- `elem_mat_nmatrices`= number of element matrices,
- `ndofs_row`= Total number of equations in the row direction,
- `ndofs_col`= Total number of equations in the column direction.

If the `buffer_pointer` field of the assembler is 0, which is the case after
that assembler was created, the buffers are resized appropriately given the
dimensions on input. Otherwise, the buffers are left completely untouched.

# Returns
- `self`: the modified assembler.


!!! note

    The buffers are initially not filled with anything meaningful.
    After the assembly, only the `(self.buffer_pointer - 1)` entries
    are meaningful numbers. Beware!
"""
function startassembly!(
    self::SysmatAssemblerSparspak,
    elem_mat_nrows,
    elem_mat_ncols,
    elem_mat_nmatrices,
    ndofs_row,
    ndofs_col;
    force_init = false
)
    # Only resize the buffers if the pointer is less than 1
    if self.ndofs_col < 1 || self.ndofs_row < 1
        self.ndofs_row = ndofs_row
        self.ndofs_col = ndofs_col
        self.problem = Problem(self.ndofs_row, self.ndofs_col, elem_mat_nrows * self.ndofs_col, zero(eltype(self.problem.rhs)))
    end
    # Leave the buffers uninitialized, unless the user requests otherwise
    if force_init

    end

    return self
end

"""
    assemble!(
        self::SysmatAssemblerSparspak,
        mat::MT,
        dofnums_row::IT,
        dofnums_col::IT,
    ) where {MT, IT}

Assemble a rectangular matrix.
"""
function assemble!(
    self::SysmatAssemblerSparspak,
    mat::MT,
    dofnums_row::IT,
    dofnums_col::IT,
) where {MT, IT}
    # Assembly of a rectangular matrix.
    # The method assembles a rectangular matrix using the two vectors of
    # equation numbers for the rows and columns.
    nrows = length(dofnums_row)
    ncolumns = length(dofnums_col)
    @assert size(mat) == (nrows, ncolumns)
    @inbounds for j in 1:ncolumns
        if 0 < dofnums_col[j] <= self.ndofs_col
            c = dofnums_col[j]
            for i in 1:nrows
                r = dofnums_row[i]
                if 0 < dofnums_row[i] <= self.ndofs_row
                    inaij!(self.problem, r, c, mat[i, j])
                end
            end
        end
    end
    return self
end

"""
    makematrix!(self::SysmatAssemblerSparspak)

Make a sparse matrix.
"""
function makematrix!(self::SysmatAssemblerSparspak{IT, FT}) where {IT, FT}
    @info "in makematrix!"
    s = SparseSolver(self.problem)
    if !self.nomatrixresult
        # @time     findorder!(s) || ErrorException("Finding Order.")
        # @time     symbolicfactor!(s) || ErrorException("Symbolic Factorization.")
        # @time     inmatrix!(s) || ErrorException("Matrix input.")
    end
    return s
end
