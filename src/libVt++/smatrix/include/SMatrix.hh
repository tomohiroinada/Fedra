#ifndef __SMATRIX_HH
#define __SMATRIX_HH
// ********************************************************************
//
// source:
//
// type:      source code
//
// created:   20. Mar 2001
//
// author:    Thorsten Glebe
//            HERA-B Collaboration
//            Max-Planck-Institut fuer Kernphysik
//            Saupfercheckweg 1
//            69117 Heidelberg
//            Germany
//            E-mail: T.Glebe@mpi-hd.mpg.de
//
// Description: A fixed size two dimensional Matrix class
//
// changes:
// 20 Mar 2001 (TG) creation
// 21 Mar 2001 (TG) added operators +=, -=, *=, /=
// 26 Mar 2001 (TG) place_in_row(), place_in_col() added
// 02 Apr 2001 (TG) non-const Array() added
// 03 Apr 2001 (TG) invert() added
// 07 Apr 2001 (TG) CTOR from SVertex (dyadic product) added
// 09 Apr 2001 (TG) CTOR from array added
// 11 Apr 2001 (TG) rows(), cols(), size() replaced by rows, cols, size
// 25 Mai 2001 (TG) row(), col() added
// 04 Sep 2001 (TG) moved inlined functions to .icc file
// 11 Jan 2002 (TG) added operator==(), operator!=()
// 14 Jan 2002 (TG) added more operator==(), operator!=(), operator>(), operator<()
//
// ********************************************************************
#include <iosfwd>

template <class T, unsigned int D> class SVector;

// expression engine
#include "Expression.hh"

/** SMatrix.
    A generic fixed size n x m Matrix class.q
    
    @memo SMatrix
    @author T. Glebe
*/
//==============================================================================
// SMatrix: column-wise storage
//==============================================================================
template <class T, unsigned int D1, unsigned int D2 = D1>
class SMatrix {
public:
  /** @name --- Typedefs --- */
  ///
  typedef T  value_type;
  
  /** @name --- Constructors --- */
  ///
  SMatrix();
  ///
  SMatrix(const SMatrix<T,D1,D2>& rhs);
  ///
  template <class A>
  SMatrix(const Expr<A,T,D1,D2>& rhs);
  /// 2nd arg: set only diagonal?
  SMatrix(const T& rhs, const bool diagonal=false);
  /// constructor via dyadic product
  SMatrix(const SVector<T,D1>& rhs);
  /// constructor via dyadic product
  template <class A>
  SMatrix(const Expr<A,T,D1>& rhs);
  /** constructor via array, triag=true: array contains only upper/lower
      triangular part of a symmetric matrix, len: length of array */
  template <class T1>
  SMatrix(const T1* a, const bool triang=false, const unsigned int len=D1*D2);

  ///
  SMatrix<T,D1,D2>& operator=(const T& rhs);
  ///
  template <class A>
  SMatrix<T,D1,D2>& operator=(const Expr<A,T,D1,D2>& rhs);

  /// return no. of matrix rows
  static const unsigned int rows = D1;
  /// return no. of matrix columns
  static const unsigned int cols = D2;
  /// return no of elements: rows*columns
  static const unsigned int size = D1*D2;

  /** @name --- Access functions --- */
  /// access the parse tree
  T apply(unsigned int i) const;
  /// return read-only pointer to internal array
  const T* Array() const;
  /// return pointer to internal array
  T* Array();

  /** @name --- Operators --- */
  /// element wise comparison
  bool operator==(const T& rhs) const;
  /// element wise comparison
  bool operator!=(const T& rhs) const;
  /// element wise comparison
  bool operator==(const SMatrix<T,D1,D2>& rhs) const;
  /// element wise comparison
  bool operator!=(const SMatrix<T,D1,D2>& rhs) const;
  /// element wise comparison
  template <class A>
  bool operator==(const Expr<A,T,D1,D2>& rhs) const;
  /// element wise comparison
  template <class A>
  bool operator!=(const Expr<A,T,D1,D2>& rhs) const;

  /// element wise comparison
  bool operator>(const T& rhs) const;
  /// element wise comparison
  bool operator<(const T& rhs) const;
  /// element wise comparison
  bool operator>(const SMatrix<T,D1,D2>& rhs) const;
  /// element wise comparison
  bool operator<(const SMatrix<T,D1,D2>& rhs) const;
  /// element wise comparison
  template <class A>
  bool operator>(const Expr<A,T,D1,D2>& rhs) const;
  /// element wise comparison
  template <class A>
  bool operator<(const Expr<A,T,D1,D2>& rhs) const;

  /// read-only access
  const T& operator()(unsigned int i, unsigned int j) const;
  /// read/write access
  T& operator()(unsigned int i, unsigned int j);
  ///
  SMatrix<T,D1,D2>& operator+=(const SMatrix<T,D1,D2>& rhs);
  ///
  template <class A>
  SMatrix<T,D1,D2>& operator+=(const Expr<A,T,D1,D2>& rhs);
  ///
  SMatrix<T,D1,D2>& operator-=(const SMatrix<T,D1,D2>& rhs);
  ///
  template <class A>
  SMatrix<T,D1,D2>& operator-=(const Expr<A,T,D1,D2>& rhs);
  ///
  SMatrix<T,D1,D2>& operator*=(const SMatrix<T,D1,D2>& rhs);
  ///
  template <class A>
  SMatrix<T,D1,D2>& operator*=(const Expr<A,T,D1,D2>& rhs);
  ///
  SMatrix<T,D1,D2>& operator/=(const SMatrix<T,D1,D2>& rhs);
  ///
  template <class A>
  SMatrix<T,D1,D2>& operator/=(const Expr<A,T,D1,D2>& rhs);

  /** @name --- Expert functions --- */
  /// invert symmetric, pos. def. Matrix via Dsinv
  bool sinvert();
  /** determinant of symmetrc, pos. def. Matrix via Dsfact. \textbf{Note:} this
      will destroy the contents of the Matrix!*/
  bool sdet(T& det);
  /// invert square Matrix via Dinv
  bool invert();
  /** determinant of square Matrix via Dfact. \textbf{Note:} this will destroy
      the contents of the Matrix! */
  bool det(T& det);
  /// place a vector in a Matrix row
  template <unsigned int D>
  SMatrix<T,D1,D2>& place_in_row(const SVector<T,D>& rhs,
				 const unsigned int row,
				 const unsigned int col);
  /// place a vector expression in a Matrix row
  template <class A, unsigned int D>
  SMatrix<T,D1,D2>& place_in_row(const Expr<A,T,D>& rhs,
				 const unsigned int row,
				 const unsigned int col);
  /// place a vector in a Matrix column
  template <unsigned int D>
  SMatrix<T,D1,D2>& place_in_col(const SVector<T,D>& rhs,
				 const unsigned int row,
				 const unsigned int col);
  /// place a vector expression in a Matrix column
  template <class A, unsigned int D>
  SMatrix<T,D1,D2>& place_in_col(const Expr<A,T,D>& rhs,
				 const unsigned int row,
				 const unsigned int col);
  /// place a matrix in this matrix
  template <unsigned int D3, unsigned int D4>
  SMatrix<T,D1,D2>& place_at(const SMatrix<T,D3,D4>& rhs,
			     const unsigned int row,
			     const unsigned int col);
  /// place a matrix expression in this matrix
  template <class A, unsigned int D3, unsigned int D4>
  SMatrix<T,D1,D2>& place_at(const Expr<A,T,D3,D4>& rhs,
			     const unsigned int row,
			     const unsigned int col);
  /// return a Matrix row as a vector
  SVector<T,D2> row(const unsigned int therow) const;
  /// return a Matrix column as a vector
  SVector<T,D1> col(const unsigned int thecol) const;
  /// used by operator<<()
  std::ostream& print(std::ostream& os) const;

private:
  T array[D1*D2];
}; // end of class SMatrix

//==============================================================================
// operator<<
//==============================================================================
template <class T, unsigned int D1, unsigned int D2>
inline std::ostream& operator<<(std::ostream& os, const SMatrix<T,D1,D2>& rhs) {
  return rhs.print(os);
}

#include "SMatrix.icc"
// include Matrix-Vector multiplication
#include "MatrixFunctions.hh"

#endif
