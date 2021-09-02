//
// vecmat3.h -  Fast 3d vector and matrix classes using expression templates
//
// Copyright (c) 2007-2013  Ramses van Zon
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.
//
// NOTES: 
//
// - This is version 2.3.1 of this header-only library. It was
//   released in May 2013, but is essentially version 2.3 from
//   February 2011, however, with an explicit open source software
//   license (MIT).
//
// - Documentation can be found in vecmat3.pdf.
//
// - The curious may find the main principles of expression templates
//   at the end of the file.
//

#ifndef _VECMAT3_
#define _VECMAT3_

// 
// Handle nonstandard c++ compiler cxx on alpha
//
#if defined(__alpha)
# include <math.h>
# include <iostream.h>
# define OSTREAM ostream
#else
# include <cmath>
# include <iostream>
# define OSTREAM std::ostream
#endif

//
// For g++ and icpc, INLINE forces inlining, even without optimization.
// In all other cases, INLINE=inline, and inlining may not occur.
// Note for xlC: 
//    In version 10, you need "-O4" to get full inlining.
//    In version 11, "-O2 -qinline=level=6" suffices.
//
#if not defined(INLINE)
# if defined(__INTEL_COMPILER)
#   define INLINE  __forceinline
# elif defined(__GNUC__)
#   define INLINE inline __attribute__((always_inline)) 
# else
#   define INLINE inline
# endif
#endif

//
// Some useful short-hand notational macros (while macros are evil,
// so is c++'s template notation, and a few abbreviations will
// alleviate our notational pain)
//
#define ENODE(A,B,C)                    typename A,int B,typename C
#define EXPRESSION_TEMPLATE             template <typename T,ENODE(X,Y,Z)>
#define EXPRESSION_TEMPLATE_PAIR        template <typename T,ENODE(A,B,C),ENODE(D,E,F)> 
#define EXPRESSION_TEMPLATE_TRIPLET     template <typename T,ENODE(A,B,C),ENODE(D,E,F),ENODE(G,H,I)>
#define EXPRESSION_TEMPLATE_MEMBER      template <ENODE(X,Y,Z)>
#define ANOTHER_EXPRESSION_TEMPLATE     template <ENODE(U,V,W)>
#define CONVERTIBLE_TEMPLATE            template <typename CONVERT>
#define CONVERTIBLE_EXPRESSION_TEMPLATE template <typename T,ENODE(X,Y,Z),typename CONVERT>
#define CONVERTIBLE_ANOTHER_EXPRESSION_TEMPLATE template <ENODE(U,V,W),typename CONVERT>
#define VECTOR           Vector<T,X,Y,Z>
#define ANOTHER_VECTOR   Vector<T,U,V,W>
#define VECTOR1          Vector<T,A,B,C>
#define VECTOR2          Vector<T,D,E,F>
#define VECTOR3          Vector<T,G,H,I>
#define MATRIX           Matrix<T,X,Y,Z>
#define ANOTHER_MATRIX   Matrix<T,U,V,W>
#define MATRIX1          Matrix<T,A,B,C>
#define MATRIX2          Matrix<T,D,E,F>
#define TT               T,Base,NoOp,Base

//
//  A dedicated namespace for all vector and matrix classes and operations.
//
namespace vecmat3 {

    //
    // Enumeration type to tag template operators
    //
    enum {
        NoOp,        // do nothing
        RowOp,       // row operation (matrix only)
        ColOp,       // column operation (matrix only)
        PlusOp,      // addition operator    
        MinusOp,     // subtraction operator
        TimesOp,     // multiplication operator
        NegativeOp,  // negation operator
        TransposeOp, // matrix transpose (matrix only)
        DyadicOp,    // dyadic of two vectors (vectors only)
        USER         // for user-defined template operations
    };
    
    //
    // A few forward declarations are needed.
    //
    class Base;  
    template <typename T,typename A=Base,int B=NoOp,typename C=Base> class Vector;
    template <typename T,typename A=Base,int B=NoOp,typename C=Base> class Matrix;
    template <typename T> class CommaOp;
    
    //
    // Default vector class
    //
    template <typename T>
    class Vector<TT> 
    {
      public:
        T x, y, z;   // vector elements
        
        //
        //  Constructors
        //
        EXPRESSION_TEMPLATE_MEMBER 
        INLINE            Vector( const VECTOR& v );
        INLINE            Vector() {}  
        INLINE explicit   Vector( const T x_, const T y_ = 0, const T z_ = 0 ); 
      
        //
        // Assignment operators
        // (inline to appease xlC compiler's issues with templated return types)
        //
        EXPRESSION_TEMPLATE_MEMBER 
        INLINE const Vector<TT>& operator= ( const VECTOR & v ) 
        {
            // No (this!=&m) clause is needed because this is a default vector,
            // and the default vector case is specialized below.
            T yValue = v.eval<1>();
            T zValue = v.eval<2>();
            x = v.eval<0>();
            y = yValue;
            z = zValue;
            return *this;
        }	
        // Specialize for a default vector
        INLINE const Vector<TT>& operator=(const Vector<TT>& v)  
        {
            if (this != &v) {
                x = v.x;
                y = v.y;
                z = v.z;
            }
            return *this;
        }
        INLINE CommaOp<T>        operator= ( const T e );// start comma operator assignment
      
        //
        //  Non-template member functions
        //  
        INLINE T          nrm()  const;  // return the norm of the 3d vector
        INLINE T          nrm2() const;  // return the squared norm       
        INLINE void       zero();        // set this vector to zero
      
        //
        //  Template evaluation (the work horse of TE)
        //
        template <int I> 
        INLINE T          eval() const;
      
        //
        //  Operators
        //
        INLINE const T&   operator() ( const int i ) const; // immutable access
        INLINE T&         operator() ( const int i );       // assignable access

        // New in 2.3: c-like square bracket support, so that v[0]=v.x, v[1]=v.y, v[2]=v.z
        INLINE const T&   operator[] ( const int i ) const;     
        INLINE T&         operator[] ( const int i );

        //
        // Template operators
        // (inline to appease xlC compiler's issues with templated return types)
        // Note: CONVERT should be a T or something that can be converted to T.
        //  
        CONVERTIBLE_TEMPLATE 
        INLINE Vector<TT>& operator*= ( const CONVERT a )
        {
            x *= a;
            y *= a;
            z *= a;
            return *this;
        }      
        CONVERTIBLE_TEMPLATE 
        INLINE Vector<TT>& operator/= ( const CONVERT a )
        {
            x /= a;
            y /= a;
            z /= a;
            return *this;
        }        
        EXPRESSION_TEMPLATE_MEMBER 
        INLINE Vector<TT>& operator+= ( const VECTOR& v )
        {
            T yValue = v.eval<1>();
            T zValue = v.eval<2>();
            x += v.eval<0>();
            y += yValue;
            z += zValue;
            return *this;
        }
        EXPRESSION_TEMPLATE_MEMBER 
        INLINE Vector<TT>& operator-= ( const VECTOR& v )
        {
            T yValue = v.eval<1>();
            T zValue = v.eval<2>();        
            x -= v.eval<0>();
            y -= yValue;
            z -= zValue;
            return *this;
        }
	
    }; // end definition vector class
    
    //
    //  The default matrix class
    //
    template <typename T> 
    class Matrix<TT> 
    {
      public:
        T xx, xy, xz;  // elements of the matrix
        T yx, yy, yz;
        T zx, zy, zz;
        
        //
        // Constructors
        //
        EXPRESSION_TEMPLATE_MEMBER 
        INLINE            Matrix( const MATRIX& m );
        INLINE            Matrix() {}
        INLINE explicit   Matrix( T a,   T b=0, T c=0,
                                  T d=0, T e=0, T f=0,
                                  T g=0, T h=0, T i=0);
        //
        // Assignment operators
        // (inline to appease xlC compiler's issues with templated return types)
        //
        EXPRESSION_TEMPLATE_MEMBER 
        INLINE const Matrix<TT>& operator= ( const MATRIX& m )
        {
            // no (this!=&m) clause needed because this is a default matrix,
            // and the default matrix case is specialized below.
            T xxValue = m.template eval<0,0>(); 
            T xyValue = m.template eval<0,1>(); 
            T xzValue = m.template eval<0,2>();
            T yxValue = m.template eval<1,0>(); 
            T yyValue = m.template eval<1,1>(); 
            T yzValue = m.template eval<1,2>();
            T zxValue = m.template eval<2,0>(); 
            T zyValue = m.template eval<2,1>();
            zz = m.template eval<2,2>();
            xx = xxValue; xy = xyValue; xz = xzValue;
            yx = yxValue; yy = yyValue; yz = yzValue;
            zx = zxValue; zy = zyValue;
            return *this;
        }
        // specialize for default matrix
        INLINE const Matrix<TT>& operator=( const Matrix<TT>& m )
        {
            if (this != &m) { 
                xx = m.xx; yx = m.yx; zx = m.zx;
                xy = m.xy; yy = m.yy; zy = m.zy;
                xz = m.xz; yz = m.yz; zz = m.zz;
            }
            return *this;
        }
        INLINE CommaOp<T>       operator=(T e);           // assign from comma expression
   
        //
        // Evaluate components
        //
        template <int I,int J> INLINE T eval() const; 

        //
        // Non-template member functions for getting individual rows or columns
        // (inline to appease xlC compiler's issues with templated return types)
        //    
        INLINE Vector<T,Matrix<TT>,RowOp,Base> row   ( const int i ) const
        {
            return Vector<T,Matrix<TT>,RowOp,Base>(*this,i);
        }        
        INLINE Vector<T,Matrix<TT>,ColOp,Base> column( const int j ) const
        {
            return Vector<T,Matrix<TT>,ColOp,Base>(*this,j);
        }
        
        // New in 2.3: c-like square bracket support, so that v[0]=v.x, v[1]=v.y, v[2]=v.z
        // inlined here, because xlC has issues with templated return types.
        INLINE Vector<TT>& operator[] ( const int i )
        {
            return *((Vector<TT>*)(&xx+3*i));
        }

        INLINE T         nrm2() const;       // 2-norm squared
        INLINE T         nrm()  const;       // norm (=square root of nrm2)
        INLINE T         tr()   const;       // trace
        INLINE T         det()  const;       // determinant
        INLINE void      zero();             // set this matrix to zero
        INLINE void      one();              // set this matrix to the identity matrix
        INLINE void      reorthogonalize();  // make this matrix is orthogonal

        //
        // Template member functions to set rows and columns
        //
        EXPRESSION_TEMPLATE_MEMBER 
        INLINE void       setRow   ( const int i, const VECTOR& v );

        EXPRESSION_TEMPLATE_MEMBER 
        INLINE void       setColumn( const int j, const VECTOR& v );

        //
        //  Non-template operators to access elements
        //
        INLINE const T&   operator() ( const int i, const int j ) const;
        INLINE T&         operator() ( const int i, const int j );

        //
        // Template operators to add, subtract, multiply and divide by.
        // CONVERT is a type T or something that can automatically be converted to T
        // (inline to appease xlC compiler's issues with templated return types)
        //
        CONVERTIBLE_TEMPLATE 
        INLINE Matrix<TT>& operator*= ( const CONVERT a )
        {
            xx *= a; xy *= a; xz *= a;
            yx *= a; yy *= a; yz *= a;
            zx *= a; zy *= a; zz *= a;
            return *this;
        }

        CONVERTIBLE_TEMPLATE 
        INLINE Matrix<TT>& operator/= ( const CONVERT a )
        {
            xx /= a; xy /= a; xz /= a;
            yx /= a; yy /= a; yz /= a;
            zx /= a; zy /= a; zz /= a;
            return *this;
        }

        EXPRESSION_TEMPLATE_MEMBER 
        INLINE Matrix<TT>& operator+= ( const MATRIX& m )
        {
            T xxValue = m.template eval<0,0>(); 
            T xyValue = m.template eval<0,1>(); 
            T xzValue = m.template eval<0,2>();
            T yxValue = m.template eval<1,0>(); 
            T yyValue = m.template eval<1,1>(); 
            T yzValue = m.template eval<1,2>();
            T zxValue = m.template eval<2,0>(); 
            T zyValue = m.template eval<2,1>(); 
            zz += m.template eval<2,2>();
            xx += xxValue; xy += xyValue; xz += xzValue;
            yx += yxValue; yy += yyValue; yz += yzValue;
            zx += zxValue; zy += zyValue;
            return *this;
        }

        EXPRESSION_TEMPLATE_MEMBER 
        INLINE Matrix<TT>& operator-= ( const MATRIX& m )
        {
            T xxValue = m.template eval<0,0>(); 
            T xyValue = m.template eval<0,1>(); 
            T xzValue = m.template eval<0,2>();
            T yxValue = m.template eval<1,0>(); 
            T yyValue = m.template eval<1,1>(); 
            T yzValue = m.template eval<1,2>();
            T zxValue = m.template eval<2,0>(); 
            T zyValue = m.template eval<2,1>(); 
            zz -= m.template eval<2,2>();
            xx -= xxValue; xy -= xyValue; xz -= xzValue;
            yx -= yxValue; yy -= yyValue; yz -= yzValue;
            zx -= zxValue; zy -= zyValue;
            return *this;
        }

        EXPRESSION_TEMPLATE_MEMBER 
        INLINE Matrix<TT>& operator*= ( const MATRIX& m )
        {  
            Matrix<TT> rhs(m);
            T x, y;
            x = xx; y = xy;
            xx *= rhs.xx; xx += y*rhs.yx+xz*rhs.zx;
            xy *= rhs.yy; xy += x*rhs.xy+xz*rhs.zy;
            xz *= rhs.zz; xz += x*rhs.xz+ y*rhs.yz;
            x = yx; y = yy;
            yx *= rhs.xx; yx += y*rhs.yx+yz*rhs.zx;
            yy *= rhs.yy; yy += x*rhs.xy+yz*rhs.zy;
            yz *= rhs.zz; yz += x*rhs.xz+ y*rhs.yz;
            x = zx; y = zy;
            zx *= rhs.xx; zx += y*rhs.yx+zz*rhs.zx;
            zy *= rhs.yy; zy += x*rhs.xy+zz*rhs.zy;
            zz *= rhs.zz; zz += x*rhs.xz+ y*rhs.yz;
            return *this;          
        }

    }; // end default matrix class definition

    //
    // Declarations of vector-matrix operations
    //
    EXPRESSION_TEMPLATE_PAIR 
    INLINE Vector<T,VECTOR1,PlusOp,VECTOR2> 
    operator+ ( const VECTOR1 & v1, 
                const VECTOR2 & v2 );

    EXPRESSION_TEMPLATE_PAIR 
    INLINE Vector<T,VECTOR1,MinusOp,VECTOR2> 
    operator- ( const VECTOR1 & v1, 
                const VECTOR2 & v2 );

    // Template expression for cross product of two vector expressions
    EXPRESSION_TEMPLATE_PAIR 
    INLINE Vector<T,VECTOR1,TimesOp,VECTOR2> 
    operator^ ( const VECTOR1 & v1, 
                const VECTOR2 & v2 );

    CONVERTIBLE_EXPRESSION_TEMPLATE 
    INLINE Vector<T,VECTOR,TimesOp,T> 
    operator/ ( const VECTOR & v,  
                CONVERT a );

    CONVERTIBLE_EXPRESSION_TEMPLATE 
    INLINE Vector<T,VECTOR,TimesOp,T> 
    operator* ( CONVERT a, 
                const VECTOR & v );

    CONVERTIBLE_EXPRESSION_TEMPLATE 
    INLINE Vector<T,VECTOR,TimesOp,T> 
    operator* ( const VECTOR& v, 
                CONVERT a );

    EXPRESSION_TEMPLATE 
    INLINE Vector<T,VECTOR,NegativeOp,Base> 
    operator- ( const VECTOR & v );

    // Inner product of two vector expressions implemented as operator|
    EXPRESSION_TEMPLATE_PAIR 
    INLINE T operator| ( const VECTOR1 & v1, 
                         const VECTOR2 & v2 );

    // Inner product of two vector expressions implemented as operator*
    EXPRESSION_TEMPLATE_PAIR 
    INLINE T operator* ( const VECTOR1 & v1, 
                         const VECTOR2 & v2);

    EXPRESSION_TEMPLATE_PAIR 
    INLINE Matrix<T,MATRIX1,PlusOp,MATRIX2> 
    operator+ ( const MATRIX1 & m1, 
                const MATRIX2 & m2 );

    EXPRESSION_TEMPLATE_PAIR 
    INLINE Matrix<T,MATRIX1,MinusOp,MATRIX2> 
    operator- ( const MATRIX1 & m1, 
                const MATRIX2 & m2 );

    EXPRESSION_TEMPLATE_PAIR 
    INLINE Matrix<T,MATRIX1,TimesOp,MATRIX2> 
    operator* ( const MATRIX1 & m1, 
                const MATRIX2 & m2);

    EXPRESSION_TEMPLATE_PAIR 
    INLINE Vector<T,MATRIX1,TimesOp,VECTOR2> 
    operator* ( const MATRIX1 & m, 
                const VECTOR2 & v );

    CONVERTIBLE_EXPRESSION_TEMPLATE 
    INLINE Matrix<T,MATRIX,TimesOp,T> 
    operator* ( CONVERT a, 
                const MATRIX & m);

    CONVERTIBLE_EXPRESSION_TEMPLATE 
    INLINE Matrix<T,MATRIX,TimesOp,T> 
    operator* ( const MATRIX & m, 
                CONVERT a );

    CONVERTIBLE_EXPRESSION_TEMPLATE 
    INLINE Matrix<T,MATRIX,TimesOp,T> 
    operator/ ( const MATRIX & m, 
                CONVERT a );

    EXPRESSION_TEMPLATE 
    INLINE Matrix<T,MATRIX,NegativeOp,Base> 
    operator- ( const MATRIX & m );

    // Matrix expression for the transpose of m
    EXPRESSION_TEMPLATE 
    INLINE Matrix<T,MATRIX,TransposeOp,Base> 
    Transpose ( const MATRIX & m );

    // Matrix expression for the dyadic matrix formed by v1 and v2
    EXPRESSION_TEMPLATE_PAIR 
    INLINE Matrix<T,VECTOR1,DyadicOp,VECTOR2> 
    Dyadic ( const VECTOR1 & v1, 
             const VECTOR2 & v2 );

    // Distance between two vectors
    EXPRESSION_TEMPLATE_PAIR 
    INLINE T dist( const VECTOR1 & v1, const VECTOR2 & v2 ); 

    // Squared distance between two vectors
    EXPRESSION_TEMPLATE_PAIR 
    INLINE T dist2( const VECTOR1 & v1, 
                    const VECTOR2 & v2 );

    // Distance between two vectors with the first one shifted by v3
    EXPRESSION_TEMPLATE_TRIPLET 
    INLINE T distwithshift( const VECTOR1 & v1, 
                            const VECTOR2 & v2, 
                            const VECTOR3 & v3 );

    // Compute the inverse of matrix m
    EXPRESSION_TEMPLATE 
    INLINE Matrix<TT> Inverse( const MATRIX & m );

    // return the rotation matrix around vector v by an angle v.nrm()
    EXPRESSION_TEMPLATE 
    INLINE Matrix<TT> Rodrigues ( const VECTOR & v );

    // output to streams for vector and matrix expressions
    EXPRESSION_TEMPLATE 
    INLINE OSTREAM& operator<< ( OSTREAM & o, const VECTOR & v );
    EXPRESSION_TEMPLATE 
    INLINE OSTREAM& operator<< ( OSTREAM & o, const MATRIX & v);

    //
    // Convenient macros:
    //
    #define crossProduct(a,b) ((a)^(b))
    #define dotProduct(a,b)   ((a)|(b))
    #define MTVmult(M,v)      (Transpose(M)*(v))

    //
    // Deprecated macros for backward compatibility:
    // (not recommended for new applications)
    //
    #define returnVector(x,y,z)                       return Vector(x,y,z)
    #define return_Vector(x,y,z)                      return Vector(x,y,z)
    #define ZERO(v)                                   v(0)
    #define Vector_ZERO(v)                            Vector v(0)
    #define Vector_INIT(v,x,y,z)                      Vector v(x,y,z)
    #define DIFF(one,two)                             one-two
    #define DIFFWITHSHIFT(one,two,shift)              one-two+shift
    #define Matrix_ZERO(m)                            Matrix m(0)
    #define Matrix_INIT(m,xx,xy,xz,yx,yy,yz,zx,zy,zz) Matrix m(xx,xy,xz,yx,yy,yz,zx,zy,zz)
    #define return_Matrix(xx,xy,xz,yx,yy,yz,zx,zy,zz) return Matrix(xx,xy,xz,yx,yy,yz,zx,zy,zz)
    #define returnMatrix(xx,xy,xz,yx,yy,yz,zx,zy,zz)  return Matrix(xx,xy,xz,yx,yy,yz,zx,zy,zz)

    /*****************************/
    /**  INLINE IMPLEMENTATION  **/
    /*****************************/

    //
    // Auxiliary functions
    //
    template <typename T> 
    INLINE T sqr(T x) 
    { 
        return x*x; 
    }

    template <typename T> 
    INLINE T absmax(const T* begin, const T* last) 
    {
        const T* p = begin;
        T max = (*p)<0?-(*p):(*p);
        for (++p; p<=last; ++p) {
            T absval = (*p)<0?-(*p):(*p);
            if (absval > max) 
            max = absval;
        }     
        return max;
    }

    //
    // Constructors
    //
    template <typename T> 
    INLINE Vector<TT>::Vector(T x_, T y_, T z_) : 
      x(x_), y(y_), z(z_) 
    {}

    template <typename T> 
    EXPRESSION_TEMPLATE_MEMBER 
    INLINE Vector<TT>::Vector(const VECTOR& v) :
      x(v.eval<0>()), y(v.eval<1>()), z(v.eval<2>())
    {}

    //
    // Member functions of the basic Vector type
    //

    // Assign zero to all elements
    template<typename T> 
    INLINE void Vector<TT>::zero() 
    {
        x = y = z = 0;
    }

    // Access to elements through (square) parentheses for assignment 
    template <typename T>
    INLINE T & Vector<TT>::operator() ( const int i ) 
    { 
        return *(&x+i);
    }

    // Passive access to elements through parentheses 
    template <typename T> 
    INLINE const T& Vector<TT>::operator() ( const int i ) const 
    { 
        return *(&x+i);
    }

    // New in 2.3: optional c-like square bracket support, so that
    // v[0]=v.x, v[1]=v.y, v[2]=v.z, and m[0][0]=m.xx ... m[2][2]=m.zz
    template <typename T> 
    INLINE T & Vector<TT>::operator[] ( const int i ) 
    { 
      return *(&x+i);
    }
    template <typename T> 
    INLINE const T & Vector<TT>::operator[] ( const int i ) const 
    { 
        return *(&x+i);
    }
    //

    // Evaluate elements (basis of the template evaluations technique)
    template <typename T> 
     template <int I> 
    INLINE T Vector<TT>::eval() const 
    { 
        switch(I) {
        case 0: return x; 
        case 1: return y; 
        case 2: return z; 
        default: return 0;
        }
    }

    // Norm squared
    template <typename T> 
    INLINE T Vector<TT>::nrm2() const 
    {
        return x*x+y*y+z*z;
    }

    // Norm
    template <typename T> 
    INLINE T Vector<TT>::nrm() const 
    {
        T max = absmax(&x,&z);
        if (max != 0)
            return max*(T)(sqrt(sqr(x/max)
                                +sqr(y/max)
                                +sqr(z/max)));
        else
            return 0;
    }

    //                                           
    // Member functions of the basic Matrix type: 
    //                                           

    //
    // Constructors
    //
    template <typename T> 
    INLINE Matrix<TT>::Matrix( T a, T b, T c,
                               T d, T e, T f,
                               T g, T h, T i ) :
      xx(a), xy(b), xz(c),
      yx(d), yy(e), yz(f),
      zx(g), zy(h), zz(i)
    {}
    template <typename T> 
    EXPRESSION_TEMPLATE_MEMBER 
    INLINE Matrix<TT>::Matrix( const MATRIX & m ) :
      xx(m.template eval<0,0>()), xy(m.template eval<0,1>()), xz(m.template eval<0,2>()),
      yx(m.template eval<1,0>()), yy(m.template eval<1,1>()), yz(m.template eval<1,2>()),
      zx(m.template eval<2,0>()), zy(m.template eval<2,1>()), zz(m.template eval<2,2>())
    {}


    // Set all elements to zero
    template <typename T> 
    INLINE void Matrix<TT>::zero() 
    {
        xx = xy = xz = yx = yy = yz = zx = zy = zz = 0;
    }

    // Turn this into an identity matrix
    template <typename T> 
    INLINE void Matrix<TT>::one() 
    {
        xx = yy = zz = 1;
        xy = xz = yx = yz = zx = zy = 0;
    }

    // Set the elements of a row equal to those of a vector
    template <typename T> 
    EXPRESSION_TEMPLATE_MEMBER 
    INLINE void Matrix<TT>::setRow( const int i, 
                                    const VECTOR & v ) 
    {
        switch(i) {
        case 0: xx = v.template eval<0>(); xy = v.template eval<1>(); xz = v.template eval<2>(); break;
        case 1: yx = v.template eval<0>(); yy = v.template eval<1>(); yz = v.template eval<2>(); break;
        case 2: zx = v.template eval<0>(); zy = v.template eval<1>(); zz = v.template eval<2>(); break;
        default: break;
        }
    }

    // Set the elements of a column equal to those of a vector
    template <typename T> 
    EXPRESSION_TEMPLATE_MEMBER 
    INLINE void Matrix<TT>::setColumn( const int j, 
                                       const VECTOR & v ) 
    {
        switch(j) {
        case 0: 
            xx = v.template eval<0>(); 
            yx = v.template eval<1>(); 
            zx = v.template eval<2>(); 
            break;
        case 1: 
            xy = v.template eval<0>(); 
            yy = v.template eval<1>(); 
            zy = v.template eval<2>(); 
            break;
        case 2: 
            xz = v.template eval<0>(); 
            yz = v.template eval<1>(); 
            zz = v.template eval<2>(); 
            break;
        default: 
            break;
        }
    }

    // Access to elements through parentheses for assignment
    template <typename T> 
    INLINE T& Matrix<TT>::operator() ( const int i, 
                                       const int j ) 
    { 
        return *(&xx+3*i+j);
    }

    // Passive:
    template <typename T> 
    INLINE const T& Matrix<TT>::operator() ( const int i, 
                                             const int j ) const 
    { 
        return *(&xx+3*i+j);
    }

    // Passive access to the elements through eval function
    template <typename T> 
    template <int I, int J> 
    INLINE T Matrix<TT>::eval() const 
    {
        switch (I) { 
        case 0: switch (J) { 
            case 0: return xx; 
            case 1: return xy; 
            case 2: return xz; 
            default: return 0; 
            } 
        case 1: switch (J) { 
            case 0: return yx; 
            case 1: return yy; 
            case 2: return yz; 
            default: return 0; 
            } 
        case 2: switch (J) { 
            case 0: return zx; 
            case 1: return zy; 
            case 2: return zz; 
            default: return 0; 
            } 
        default: return 0;
        }
    }

    // Trace
    template <typename T> 
    INLINE T Matrix<TT>::tr() const 
    { 
        return xx + yy + zz; 
    }

    // Determinant
    template <typename T> 
    INLINE T Matrix<TT>::det() const 
    {
        return xx*(yy*zz-yz*zy)+xy*(yz*zx-yx*zz)+xz*(yx*zy-yy*zx);
    }

    // Norm squared
    template <typename T> 
    INLINE T Matrix<TT>::nrm2() const 
    { 
        return sqr(xx) + sqr(xy) + sqr(xz) + sqr(yx) + sqr(yy) 
             + sqr(yz) + sqr(zx) + sqr(zy) + sqr(zz); 
    }

    // Norm of the  matrix
    template <typename T> 
    INLINE T Matrix<TT>::nrm() const 
    {
        T max = absmax(&xx,&zz);
        if (max!= 0) 
            return max*(T)(sqrt(sqr(xx/max)+sqr(xy/max)+sqr(xz/max)
                               +sqr(yx/max)+sqr(yy/max)+sqr(yz/max)
                               +sqr(zx/max)+sqr(zy/max)+sqr(zz/max)
                                ));
        else 
            return 0;
    } 

    //
    // Member functions of Vector and Matrix:
    //

    // Define generic members functions of Vector and Matrix

    // To access the elements of a vector expression:
    #define VECPARENTHESES					\
        INLINE T operator()(int i) const {			\
            switch(i){						\
            case 0: return eval<0>();                           \
            case 1: return eval<1>();                           \
            case 2: return eval<2>();                           \
            default: return 0;                                  \
            }							\
        }                                                       \
        INLINE T operator[](int i) const {			\
            switch(i){						\
            case 0: return eval<0>();                           \
            case 1: return eval<1>();                           \
            case 2: return eval<2>();                           \
            default: return 0;                                  \
            }							\
        }
    
    // To access the elements of a matrix expression
    #define MATPARENTHESES                                      \
        INLINE T operator()(int i, int j) const {               \
            switch(i){                                          \
            case 0: switch(j){                                  \
                case 0: return eval<0,0>();                     \
                case 1: return eval<0,1>();                     \
                case 2: return eval<0,2>();                     \
                default: return 0;                              \
                }                                               \
            case 1: switch(j){                                  \
                case 0: return eval<1,0>();                     \
                case 1: return eval<1,1>();			\
                case 2: return eval<1,2>();			\
                default: return 0;                              \
                }                                               \
            case 2: switch(j){                                  \
                case 0: return eval<2,0>();			\
                case 1: return eval<2,1>();			\
                case 2: return eval<2,2>();			\
                default: return 0;				\
                }                                               \
            default:                                            \
                return 0;					\
            }							\
        }

    // To get the norm of a vector
    #define VECNRM						\
        INLINE T nrm() const {					\
            T x[3] = {eval<0>(), eval<1>(), eval<2>()};		\
            T max = absmax(x,x+2);				\
            if (max!= 0)                                        \
                return max*(T)(sqrt(sqr(x[0]/max)               \
                                    +sqr(x[1]/max)              \
                                    +sqr(x[2]/max)));		\
            else                                                \
                return 0;                                       \
        }
   
    // To get the norm of a matrix
    #define MATNRM						\
        INLINE T nrm() const {                                  \
            T x[9] = {eval<0,0>(), eval<0,1>(), eval<0,2>(),	\
                      eval<1,0>(), eval<1,1>(), eval<1,2>(),	\
                      eval<2,0>(), eval<2,1>(), eval<2,2>()};	\
            T max = absmax(x,x+8);                              \
            if (max!= 0)                                        \
                return max*(T)(sqrt(sqr(x[0]/max)               \
                                    +sqr(x[1]/max)              \
                                    +sqr(x[2]/max)              \
                                    +sqr(x[3]/max)              \
                                    +sqr(x[4]/max)              \
                                    +sqr(x[5]/max)              \
                                    +sqr(x[6]/max)              \
                                    +sqr(x[7]/max)              \
                                    +sqr(x[8]/max)));		\
            else                                                \
                return 0;                                       \
        } 

    // To get the norm squared of a vector expression
    #define VECNRM2                                             \
        INLINE T nrm2() const {                                 \
            return sqr(eval<0>())+sqr(eval<1>())+sqr(eval<2>());\
        }

    // To get the norm squared of a matrix expression
    #define MATNRM2                                                   \
        INLINE T nrm2() const {   	                              \
            return sqr(eval<0,0>())+sqr(eval<0,1>())+sqr(eval<0,2>()) \
                  +sqr(eval<1,0>())+sqr(eval<1,1>())+sqr(eval<1,2>()) \
                  +sqr(eval<2,0>())+sqr(eval<2,1>())+sqr(eval<2,2>());\
        }
    
    // To get the trace of a matrix expression
    #define MATTR                                               \
        INLINE T tr() const {					\
            return eval<0,0>() + eval<1,1>() + eval<2,2>();	\
        }

    // To get the determinant of a matrix expression
    #define MATDET                                              \
        INLINE T det() const {                                  \
            T xx=eval<0,0>(), xy=eval<0,1>(), xz=eval<0,2>();	\
            T yx=eval<1,0>(), yy=eval<1,1>(), yz=eval<1,2>();	\
            T zx=eval<2,0>(), zy=eval<2,1>(), zz=eval<2,2>();	\
            return                                              \
            xx*(yy*zz-yz*zy)+xy*(yz*zx-yx*zz)+xz*(yx*zy-yy*zx);	\
        }
    
    // To get the i-th row of a matrix expression
    #define MATROW                                              \
        INLINE Vector<T,CLASS,RowOp,Base> row(int i) const {	\
            return Vector<T,CLASS,RowOp,Base>(*this, i);        \
        }

    // To get the j-th column of a matrix expression
    #define MATCOLUMN                                           \
        INLINE Vector<T,CLASS,ColOp,Base> column(int j) const { \
            return Vector<T,CLASS,ColOp,Base>(*this, j);        \
        }

    #define VECDEFS VECNRM2 VECNRM VECPARENTHESES
    #define MATDEFS MATNRM2 MATNRM MATPARENTHESES MATDET \
                    MATTR MATROW MATCOLUMN
 
    // Expression template class for '+'  between two Vector expressions

    #define CLASS Vector<T,VECTOR1,PlusOp,VECTOR2> 
    EXPRESSION_TEMPLATE_PAIR
    class CLASS
    {
      public:
        VECDEFS
        template <int I> INLINE T eval() const;
        INLINE Vector ( const VECTOR1 & left, 
                        const VECTOR2 & right ) : 
          l(&left), 
          r(&right)
        {}
      private:
        const VECTOR1* l;  
        const VECTOR2* r;
    };

    EXPRESSION_TEMPLATE_PAIR 
    template <int I> 
    INLINE T CLASS::eval() const 
    { 
        return l->eval<I>() + r->eval<I>(); 
    }

    #undef CLASS

    // Expression template class for '-'  between two Vector expressions:

    #define CLASS Vector<T,VECTOR1,MinusOp,VECTOR2> 
    EXPRESSION_TEMPLATE_PAIR
    class CLASS
    {
      public:
        VECDEFS
        template <int I> INLINE T eval() const;
        INLINE Vector(const VECTOR1& left, const VECTOR2& right) : l(&left), r(&right) {}
      private:
        const VECTOR1* l;
        const VECTOR2* r;
    };

    EXPRESSION_TEMPLATE_PAIR 
    template <int I> 
    INLINE T CLASS::eval() const 
    {
        return l->eval<I>() - r->eval<I>(); 
    }

    #undef CLASS
  
    // Expression template class for '^' between two Vector expressions

    #define CLASS Vector<T,VECTOR1,TimesOp,VECTOR2> 
    EXPRESSION_TEMPLATE_PAIR
    class CLASS
    {
      public:
       VECDEFS
       template <int I> INLINE T eval() const;
       INLINE Vector(const VECTOR1& left, const VECTOR2& right) : l(&left), r(&right) {}
      private:
       const VECTOR1* l;  // pointers to the sub-expressions
       const VECTOR2* r;
    };

    EXPRESSION_TEMPLATE_PAIR 
    template <int I> 
    INLINE T CLASS::eval() const  
    {
       switch (I) {
          case 0:
             return l->eval<1>() * r->eval<2>() - l->eval<2>() * r->eval<1>();
             break;
          case 1:
             return l->eval<2>() * r->eval<0>() - l->eval<0>() * r->eval<2>();
             break;
          case 2:
             return l->eval<0>() * r->eval<1>() - l->eval<1>() * r->eval<0>();
             break;
          default: 
             return 0;
             break;
       }
    }

    #undef CLASS

    // Expression template class for '*' between a Vector and a T

    #define CLASS Vector<T,VECTOR,TimesOp,T> 
    EXPRESSION_TEMPLATE
    class CLASS
    {
      public:
        VECDEFS
        template <int I> INLINE T eval() const;
        CONVERTIBLE_TEMPLATE 
        INLINE Vector& operator*(CONVERT c)
        // optimize a second multiplication with T
        {
            r *= c;
            return *this;
        }
        CONVERTIBLE_TEMPLATE 
        INLINE Vector& operator/(CONVERT c)
        // optimize a subsequent division by T
        {
            r /= c;
            return *this;
        }
        INLINE Vector(const VECTOR& left, T right) : l(&left), r(right) {}
      private:
        const VECTOR* l; 
        T r;
        // befriend operators that optimize further T multiplication
        CONVERTIBLE_ANOTHER_EXPRESSION_TEMPLATE  
        friend Vector<T,ANOTHER_VECTOR,TimesOp,T>& 
            operator*(CONVERT b, Vector<T,ANOTHER_VECTOR,TimesOp,T>& a);
    };

    EXPRESSION_TEMPLATE 
    template <int I> INLINE T CLASS::eval() const 
    {
        return l->eval<I>() * r; 
    }

    #undef CLASS

    // Expression template class for unary '-' acting on a Vector

    #define CLASS Vector<T,VECTOR,NegativeOp,Base> 

    EXPRESSION_TEMPLATE
    class CLASS
    {
      public:      
        VECDEFS
        template <int I> INLINE T eval() const;
        INLINE Vector(const VECTOR& left) : l(&left) {}
      private:
        const VECTOR* l;
    };

    EXPRESSION_TEMPLATE 
    template <int I> 
    INLINE T CLASS::eval() const 
    {
        return - l->eval<I>();
    }

    #undef CLASS

    // Expression template class for '+' between two Matrix expressions:

    #define CLASS Matrix<T,MATRIX1,PlusOp,MATRIX2>

    EXPRESSION_TEMPLATE_PAIR 
    class CLASS
    {
      public:
        MATDEFS
        template <int I, int J> INLINE T eval() const;
        INLINE Matrix(const MATRIX1& left, const MATRIX2& right) : l(&left), r(&right) {}
    private:
       const MATRIX1* l;
       const MATRIX2* r;
    };

    EXPRESSION_TEMPLATE_PAIR 
    template <int I,int J> 
    INLINE T CLASS::eval() const 
    { 
       return l->eval<I,J>() + r->eval<I,J>(); 
    }

    #undef CLASS

    // Expression template class for '-' between two Matrix expressions

    #define CLASS Matrix<T,MATRIX1,MinusOp,MATRIX2> 

    EXPRESSION_TEMPLATE_PAIR 
    class CLASS
    {
    public:
       MATDEFS
       template <int I, int J> INLINE T eval() const;
       INLINE Matrix(const MATRIX1& left, const MATRIX2& right): l (&left), r(&right) {}
    private:
       const MATRIX1* l; // pointers to the sub-expressions
       const MATRIX2* r;
    };

    EXPRESSION_TEMPLATE_PAIR 
    template <int I,int J> 
    INLINE T CLASS::eval() const 
    { 
       return l->eval<I,J>() - r->eval<I,J>(); 
    }

    #undef CLASS

    // Expression template class for '*' between a Matrix and a T

    #define CLASS Matrix<T,MATRIX,TimesOp,T>

    EXPRESSION_TEMPLATE
    class CLASS
    {
    public:
       MATDEFS
       INLINE Matrix(const MATRIX& left, T right) : l(&left), r(right) {}
       template <int I, int J> INLINE T eval() const;
       CONVERTIBLE_TEMPLATE INLINE Matrix& operator*(CONVERT c);
       // further multiplication with a constant
       CONVERTIBLE_TEMPLATE INLINE Matrix& operator/(CONVERT c);
       // further division by a constant
    private:
       const MATRIX* l;
       T r;
       // be-friend operators that optimize further T multiplication
       CONVERTIBLE_ANOTHER_EXPRESSION_TEMPLATE 
       friend Matrix<T,ANOTHER_MATRIX,TimesOp,T>& operator*(CONVERT b, Matrix<T,ANOTHER_MATRIX,TimesOp,T>& a);
    };

    EXPRESSION_TEMPLATE 
    template <int I,int J> 
    INLINE T CLASS::eval() const 
    { 
       return l->eval<I,J>() * r; 
    }

    EXPRESSION_TEMPLATE 
    CONVERTIBLE_TEMPLATE 
    INLINE CLASS& CLASS::operator*(CONVERT c) 
    {
       r *= c;
       return *this;
    }

    EXPRESSION_TEMPLATE 
    CONVERTIBLE_TEMPLATE 
    INLINE CLASS& CLASS::operator/(CONVERT c) 
    {
       r /= c;
       return *this;
    }

    #undef CLASS

    // Expression template class for '*' between two Matrix expressions

    #define CLASS Matrix<T,MATRIX1,TimesOp,MATRIX2>
    EXPRESSION_TEMPLATE_PAIR 
    class CLASS
    {
      public:
        MATDEFS
        template <int I, int J> INLINE T eval() const;
        INLINE Matrix(const MATRIX1& left, const MATRIX2& right) : l (&left), r(&right) {}
      private:       
        const MATRIX1* l;
        const MATRIX2* r;
    };

    EXPRESSION_TEMPLATE_PAIR 
    template <int I, int J> 
    INLINE T CLASS::eval() const 
    {
       return l->eval<I,0>() * r->eval<0,J>()
          + l->eval<I,1>() * r->eval<1,J>()
          + l->eval<I,2>() * r->eval<2,J>();
    }

    #undef CLASS

    // Expression template class for unary '-' acting on a Matrix expression 

    #define CLASS Matrix<T,MATRIX,NegativeOp,Base>

    EXPRESSION_TEMPLATE
    class CLASS
    {
      public:
        MATDEFS
        template <int I, int J> INLINE  T eval() const;
        INLINE Matrix(const MATRIX& right) : r(&right) {}
      private:
        const MATRIX* r;
    };

    EXPRESSION_TEMPLATE 
    template <int I, int J> 
    INLINE T CLASS::eval() const 
    {
       return - r->eval<I,J>();
    }

    #undef CLASS

    // Expression template class for transpose function:

    #define CLASS Matrix<T,MATRIX,TransposeOp,Base>

    EXPRESSION_TEMPLATE
    class CLASS
    {
    public:
       MATDEFS
       template <int I, int J> INLINE T eval() const;
       INLINE Matrix (const MATRIX& left) : l(&left) {}
    private:
       const MATRIX* l;
    };

    EXPRESSION_TEMPLATE 
    template <int I, int J> 
    INLINE T CLASS::eval() const 
    {
       return l->eval<J,I>();
    }

    #undef CLASS

    // Expression template class for Dyadic operation between two Vectors

    #define CLASS Matrix<T,VECTOR1,DyadicOp,VECTOR2>

    EXPRESSION_TEMPLATE_PAIR
    class CLASS
    {
      public:
        MATDEFS
        template <int I, int J> INLINE T eval() const;
        INLINE Matrix(const VECTOR1& left, const VECTOR2& right) :l (&left), r(&right) {}
      private:
        const VECTOR1* l;
        const VECTOR2* r;
    };

    EXPRESSION_TEMPLATE_PAIR 
    template <int I, int J> 
    INLINE T CLASS::eval() const 
    {
       return l->eval<I>() * r->eval<J>();
    }

    #undef CLASS

    // Expression template class for row operation acting on a Matrix expression 

    #define CLASS Vector<T,MATRIX,RowOp,Base> 

    EXPRESSION_TEMPLATE
    class CLASS
    {
      public:
        VECDEFS
        template <int J> INLINE T eval() const;
        INLINE Vector(const MATRIX& left, int _i) : l(&left), i(_i) {}
      private:
        const MATRIX* l;
        const int  i;
    };

    EXPRESSION_TEMPLATE 
    template <int J> 
    INLINE T CLASS::eval() const 
    {
       switch(i) {
          case 0: return l->eval<0,J>();
          case 1: return l->eval<1,J>();
          case 2: return l->eval<2,J>();
          default: return 0;
       }
    }

    #undef CLASS

    // Expression template for the column op acting on a Matrix expression:

    #define CLASS Vector<T,MATRIX,ColOp,Base> 

    EXPRESSION_TEMPLATE
    class CLASS
    {
      public:       
        VECDEFS
        template <int I> INLINE T eval() const;
        INLINE Vector(const MATRIX& left, int _j) : l(&left), j(_j) {}
      private:
        const MATRIX* l;
        const int  j;
    };

    EXPRESSION_TEMPLATE 
    template <int I> 
    INLINE T CLASS::eval() const 
    {
        switch(j) {
        case 0: return l->eval<I,0>();
        case 1: return l->eval<I,1>();
        case 2: return l->eval<I,2>();
        default: return 0;
        }
    }

    #undef CLASS

    // Expression template for '*' between a Matrix and a Vector 

    #define CLASS Vector<T,MATRIX1,TimesOp,VECTOR2> 

    EXPRESSION_TEMPLATE_PAIR
    class CLASS
    {
      public:
        VECDEFS
        template <int I> INLINE T eval() const;
        INLINE Vector(const MATRIX1& left, const VECTOR2& right) : l (&left), r(&right) {}
      private:
        const MATRIX1* l;
        const VECTOR2* r;
    };

    EXPRESSION_TEMPLATE_PAIR 
    template <int I> 
    INLINE T CLASS::eval() const 
    {
       return l->eval<I,0>() * r->eval<0>()
          + l->eval<I,1>() * r->eval<1>()
          + l->eval<I,2>() * r->eval<2>();
    }

    #undef CLASS

    //
    // Helper class to implemement a comma list assignment
    //
    template <typename T>
    class CommaOp 
    {
      public:
        INLINE CommaOp& operator,(T e)  // fill a spot and go to the next
        { 
            // make sure the vector or matrix is not completely filled yet
            if (ptr <= end) 
                *ptr++ = e;  
            return *this; 
        }

        INLINE ~CommaOp()
        // fill whatever elements remain with zeros
        { 
            while (ptr <= end) 
                *ptr++ = 0;
        }

      private:
        T* ptr;       // points to the next element to be filled 
        T* const end; // points to the last element 

        INLINE CommaOp(T& e1, T& e2) :
         ptr(&e1), 
         end(&e2) 
        {}

        // Since this constructor is private, the Comma class can
        // only be used by the friend member functions operator= of
        // Vector<TT> and Matrix<TT>
        friend CommaOp<T> Vector<TT>::operator= ( T e );
        friend CommaOp<T> Matrix<TT>::operator= ( T e );
    };

    // Assignment to Vector<TT> triggers a comma separated initializer list
    template <typename T> 
    INLINE CommaOp<T> Vector<TT>::operator= ( T e ) 
    {
        x = e;
        return CommaOp<T>(y, z);
    }

    // Assignment to Matrix<TT> triggers a comma separated initializer list
    template <typename T> 
    INLINE CommaOp<T> Matrix<TT>::operator= ( T  e) 
    {
        xx = e;
        return CommaOp<T>(xy, zz);
    }

    //
    // Definitions of the operators
    //

    // Vector + Vector
    EXPRESSION_TEMPLATE_PAIR 
    INLINE Vector<T,VECTOR1,PlusOp,VECTOR2> 
    operator+ ( const VECTOR1 & v1, 
                const VECTOR2 & v2 ) 
    { 
        return Vector<T,VECTOR1,PlusOp,VECTOR2>(v1, v2);
    }

    // Vector - Vector
    EXPRESSION_TEMPLATE_PAIR 
    INLINE Vector<T,VECTOR1,MinusOp,VECTOR2> 
    operator- ( const VECTOR1 & v1, 
                const VECTOR2 & v2 ) 
    { 
        return Vector<T,VECTOR1,MinusOp,VECTOR2>(v1, v2);
    }

    // Vector ^ Vector (cross product)
    EXPRESSION_TEMPLATE_PAIR 
    INLINE Vector<T,VECTOR1,TimesOp,VECTOR2> 
    operator^ ( const VECTOR1 & v1, 
                const VECTOR2 & v2 ) 
    { 
        return Vector<T,VECTOR1,TimesOp,VECTOR2>(v1, v2);
    }

    // T * Vector
    CONVERTIBLE_EXPRESSION_TEMPLATE 
    INLINE Vector<T,VECTOR,TimesOp,T> 
    operator* ( CONVERT a, 
                const VECTOR & v ) 
    { 
        return Vector<T,VECTOR,TimesOp,T>(v, a);
    }

    // Vector * T
    CONVERTIBLE_EXPRESSION_TEMPLATE 
    INLINE Vector<T,VECTOR,TimesOp,T> 
    operator* ( const VECTOR & v, 
                CONVERT a ) 
    { 
        return Vector<T,VECTOR,TimesOp,T>(v, a);
    }

    // Vector / T
    CONVERTIBLE_EXPRESSION_TEMPLATE 
    INLINE Vector<T,VECTOR,TimesOp,T> 
    operator/ ( const VECTOR & v, 
                CONVERT a ) 
    { 
        return Vector<T,VECTOR,TimesOp,T>(v, 1/a);
    }

    // Vector * T * T
    CONVERTIBLE_EXPRESSION_TEMPLATE 
    INLINE Vector<T,VECTOR,TimesOp,T>& 
    operator* ( CONVERT a, 
                Vector<T,VECTOR,TimesOp,T> & v ) 
    { 
        v.r *= a; 
        return v; 
    }

    // Vector * T / T
    CONVERTIBLE_EXPRESSION_TEMPLATE 
    INLINE Vector<T,VECTOR,TimesOp,T>& 
    operator/ ( Vector<T,VECTOR,TimesOp,T> & v, 
                CONVERT a ) 
    { 
        v.r /= a; 
        return v; 
    }

    // - Vector
    EXPRESSION_TEMPLATE 
    INLINE Vector<T,VECTOR,NegativeOp,Base> 
    operator- ( const VECTOR & v ) 
    { 
        return Vector<T,VECTOR,NegativeOp,Base>(v);
    }

    // Vector | Vector (inner product)
    EXPRESSION_TEMPLATE_PAIR 
    INLINE T operator|( const VECTOR1 & v1,
                        const VECTOR2 & v2 ) 
    {
        return v1.template eval<0>()*v2.template eval<0>() 
             + v1.template eval<1>()*v2.template eval<1>() 
             + v1.template eval<2>()*v2.template eval<2>();
    }

    // Vector * Vector (inner product)
    EXPRESSION_TEMPLATE_PAIR 
    INLINE T 
    operator*( const VECTOR1 & v1, 
               const VECTOR2 & v2 ) 
    {
        return v1.template eval<0>()*v2.template eval<0>() 
             + v1.template eval<1>()*v2.template eval<1>() 
             + v1.template eval<2>()*v2.template eval<2>();
    }

    // Matrix + Matrix
    EXPRESSION_TEMPLATE_PAIR 
    INLINE Matrix<T,MATRIX1,PlusOp,MATRIX2> 
    operator+ ( const MATRIX1 & v1, 
                const MATRIX2 & v2 ) 
    { 
        return Matrix<T,MATRIX1,PlusOp,MATRIX2>(v1, v2);
    }

    // Matrix - Matrix
    EXPRESSION_TEMPLATE_PAIR 
    INLINE Matrix<T,MATRIX1,MinusOp,MATRIX2>
    operator-( const MATRIX1 & v1, 
               const MATRIX2 & v2 ) 
    { 
        return Matrix<T,MATRIX1,MinusOp,MATRIX2>(v1, v2);
    }

    // T * Matrix
    CONVERTIBLE_EXPRESSION_TEMPLATE 
    INLINE Matrix<T,MATRIX,TimesOp,T> 
    operator* ( CONVERT a, const MATRIX & m ) 
    { 
        return Matrix<T,MATRIX,TimesOp,T>(m, a);
    }

    // Matrix * T
    CONVERTIBLE_EXPRESSION_TEMPLATE 
    INLINE Matrix<T,MATRIX,TimesOp,T> 
    operator* ( const MATRIX & m, CONVERT a ) 
    { 
        return Matrix<T,MATRIX,TimesOp,T>(m, a);
    }

    // Matrix / T
    CONVERTIBLE_EXPRESSION_TEMPLATE 
    INLINE Matrix<T,MATRIX,TimesOp,T> 
    operator/ ( const MATRIX & m, 
                CONVERT a ) 
    { 
        return Matrix<T,MATRIX,TimesOp,T>(m, 1.0/a);
    }

    // Matrix * T * T
    CONVERTIBLE_EXPRESSION_TEMPLATE 
    INLINE Matrix<T,MATRIX,TimesOp,T>& 
    operator* ( CONVERT a, 
                Matrix<T,MATRIX,TimesOp,T> & m ) 
    { 
        m.r *= a; 
        return m; 
    }

    // Matrix * T / T
    CONVERTIBLE_EXPRESSION_TEMPLATE 
    INLINE Matrix<T,MATRIX,TimesOp,T>& 
    operator/ ( Matrix<T,MATRIX,TimesOp,T> & m,
                CONVERT a ) 
    { 
       m.r /= a; 
       return m; 
    }

    // Matrix * Matrix
    EXPRESSION_TEMPLATE_PAIR 
    INLINE Matrix<T,MATRIX1,TimesOp,MATRIX2> 
    operator* ( const MATRIX1 & a, const MATRIX2 & b ) 
    { 
        return Matrix<T,MATRIX1,TimesOp,MATRIX2>(a, b);
    }

    // -Matrix 
    EXPRESSION_TEMPLATE 
    INLINE Matrix<T,MATRIX,NegativeOp,Base> 
    operator- ( const MATRIX & m ) 
    { 
        return Matrix<T,MATRIX,NegativeOp,Base>(m);
    }

    // Matrix * Vector
    EXPRESSION_TEMPLATE_PAIR 
    INLINE Vector<T,MATRIX1,TimesOp,VECTOR2> 
    operator* ( const MATRIX1 & m, 
                const VECTOR2 & v ) 
    { 
        return Vector<T,MATRIX1,TimesOp,VECTOR2> (m, v);
    }

    // Reorthogonalize rows using Gramm-Schmidt orthogonalization of the rows
    template <typename T> 
    INLINE void Matrix<TT>::reorthogonalize() 
    {
        T z;
        int num = 10;
        while ( fabs(det()-1)> 1E-16 and --num ) {
            z = 1/row(0).nrm();   
            xx *= z;    
            xy *= z;    
            xz *= z;    
            z = xx*yx + xy*yy + xz*yz;          
            yx -= z*xx; 
            yy -= z*xy; 
            yz -= z*xz; 
            z = 1/row(1).nrm();                
            yx *= z;  
            yy *= z;  
            yz *= z;        
            zx = xy*yz - xz*yy;                 
            zy = xz*yx - xx*yz;
            zz = xx*yy - xy*yx;
            z = 1/row(2).nrm();                
            zx *= z;  
            zy *= z;  
            zz *= z;        
        }
    }

    //
    // Non-member functions
    //

    // Return |a-b|
    EXPRESSION_TEMPLATE_PAIR 
    INLINE T dist( const VECTOR1 & v1, 
                   const VECTOR2 & v2 ) 
    {
        T x[3] = {v1.template eval<0>() - v2.template eval<0>(),
                  v1.template eval<1>() - v2.template eval<1>(),
                  v1.template eval<2>() - v2.template eval<2>()};
        T max = absmax(x,x+2);  
        if (max != 0) {
            x[0] /= max;
            x[1] /= max;
            x[2] /= max;
            return max*sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
        }
        else 
            return 0;
    }

    // Return (a-b)|(a-b)
    EXPRESSION_TEMPLATE_PAIR 
    INLINE T dist2( const VECTOR1 & v1, 
                    const VECTOR2 & v2 ) 
    {
        T  d = v1.template eval<0>() - v2.template eval<0>();
        T  c = d * d;
        d = v1.template eval<1>() - v2.template eval<1>();
        c += d * d;
        d = v1.template eval<2>() - v2.template eval<2>();
        return c + d * d;
    }

    // Return |a+s-b|
    EXPRESSION_TEMPLATE_TRIPLET 
    INLINE T distwithshift( const VECTOR1 & v1, 
                            const VECTOR2 & v2, 
                            const VECTOR3 & v3 ) 
    {
        T  x = fabs(v3.template eval<0>()+v1.template eval<0>()-v2.template eval<0>());
        T  y = fabs(v3.template eval<1>()+v1.template eval<1>()-v2.template eval<1>());
        T  z = fabs(v3.template eval<2>()+v1.template eval<2>()-v2.template eval<2>());
        T  max = x;
        if (y>max) max = y;
        if (z>max) max = z;   
        if (max != 0) {
            x /= max;
            y /= max;
            z /= max;
            return max*sqrt(x*x+y*y+z*z);
        }
        else 
            return 0;
    }
    
    // Inverse of a matrix
    EXPRESSION_TEMPLATE 
    INLINE Matrix<TT> 
    Inverse(const MATRIX& me) 
    {
        Matrix<TT> m = me;
        T s = 1/m.det();
        return Matrix<TT>(
            (m.yy*m.zz-m.yz*m.zy)*s, -(m.xy*m.zz-m.xz*m.zy)*s,  (m.xy*m.yz-m.xz*m.yy)*s,
            -(m.yx*m.zz-m.yz*m.zx)*s,  (m.xx*m.zz-m.xz*m.zx)*s, -(m.xx*m.yz-m.xz*m.yx)*s,
             (m.yx*m.zy-m.yy*m.zx)*s, -(m.xx*m.zy-m.xy*m.zx)*s,  (m.xx*m.yy-m.xy*m.yx)*s);
    }

    // Build rotation matrix using the Rodrigues formula
    EXPRESSION_TEMPLATE 
    INLINE Matrix<TT> 
    Rodrigues(const VECTOR& ve) 
    {
        Vector<TT> v = ve;
        T theta = v.nrm();
        if (theta != 0) {
            T s, c;
            #ifdef SINCOS
            SINCOS(theta, &s, &c);
            #else
            s = sin(theta);
            c = cos(theta);
            #endif
            T inrm = 1/theta;
            T wx = v.x*inrm;
            T wy = v.y*inrm;
            T wz = v.z*inrm;
            T oneminusc = 1-c;
            T wxwy1mc = wx*wy*oneminusc;
            T wxwz1mc = wx*wz*oneminusc;
            T wywz1mc = wy*wz*oneminusc;
            T wxs = wx*s;
            T wys = wy*s;
            T wzs = wz*s;
            return Matrix<TT>(c+wx*wx*oneminusc, wxwy1mc-wzs,       wxwz1mc+wys,
                              wxwy1mc+wzs,       c+wy*wy*oneminusc, wywz1mc-wxs,
                              wxwz1mc-wys,       wywz1mc+wxs,       c+wz*wz*oneminusc);
        } else {
            return Matrix<TT>(1,0,0,0,1,0,0,0,1);
        }
    }

    // Transpose matrix
    EXPRESSION_TEMPLATE 
    INLINE Matrix<T,MATRIX,TransposeOp,Base> 
    Transpose( const MATRIX & m ) 
    { 
        return Matrix<T,MATRIX,TransposeOp,Base>(m);
    }

    // Dyadic product of two vectors
    EXPRESSION_TEMPLATE_PAIR 
    INLINE Matrix<T,VECTOR1,DyadicOp,VECTOR2> 
    Dyadic( const VECTOR1 & a, 
            const VECTOR2 & b ) 
    {
        return Matrix<T,VECTOR1,DyadicOp,VECTOR2> (a,b);
    }

    //
    // Output
    //

    // Vectors
    EXPRESSION_TEMPLATE 
    INLINE OSTREAM& operator<< ( OSTREAM & o, 
                                 const VECTOR & v) 
    {
        return o << v.template eval<0>() << " " 
                 << v.template eval<1>() << " " 
                 << v.template eval<2>();
    }

    // Matrices
    EXPRESSION_TEMPLATE 
    INLINE OSTREAM & operator<< ( OSTREAM & o, 
                           const MATRIX & m ) 
    {
        return o << "\n"<< m.template eval<0,0>() 
                 << " " << m.template eval<0,1>() 
                 << " " << m.template eval<0,2>()
                 << "\n"<< m.template eval<1,0>() 
                 << " " << m.template eval<1,1>() 
                 << " " << m.template eval<1,2>()
                 << "\n"<< m.template eval<2,0>() 
                 << " " << m.template eval<2,1>() 
                 << " " << m.template eval<2,2>() 
                 << "\n";
    }

}  // end namespace vecmat3

#undef VECNRM
#undef VECNRM2
#undef VECPARENTHESES
#undef VECDEFS
#undef MATNRM
#undef MATNRM2
#undef MATPARENTHESES
#undef MATTR
#undef MATDET
#undef MATROW
#undef MATCOLUMN
#undef MATDEFS
#undef OSTREAM
#undef EXPRESSION_TEMPLATE     
#undef ANOTHER_EXPRESSION_TEMPLATE 
#undef EXPRESSION_TEMPLATE_PAIR   
#undef EXPRESSION_TEMPLATE_TRIPLET
#undef EXPRESSION_TEMPLATE_MEMBER
#undef CONVERTIBLE_TEMPLATE     
#undef CONVERTIBLE_EXPRESSION_TEMPLATE     
#undef CONVERTIBLE_ANOTHER_EXPRESSION_TEMPLATE 
#undef VECTOR 
#undef VECTOR1
#undef VECTOR2
#undef VECTOR3
#undef MATRIX 
#undef MATRIX1
#undef MATRIX2
#undef TT
#undef ENODE

#ifndef NOVECMAT3DEF

// Backwards compatible type definitions of Vector and Matrix
// To not use these, state:
// #define NOVECMAT3DEF
// #include "vecmat3.h"
#ifndef DOUBLE
// default type of elements is double
#define DOUBLE double
#endif
typedef vecmat3::Vector<DOUBLE> Vector; 
typedef vecmat3::Matrix<DOUBLE> Matrix;
#endif

#endif

//
//  On the implementations of the template expressions:
// 
//  To have e.g. 'c = a+b' be unrolled to 
//  'c.x = a.x+b.x; c.y = a.y+b.y; c.z = a.z+b.z;',
//  an expression template class is defined with the following properties:
// 
//  - It contains pointers to the vectors 'a' and 'b'.
// 
//  - The constructor of this class sets up these pointers
// 
//  - It has an int template parameter indicating the type of operation.
//    The possible operations are listed in the unnamed 'enum' below.
// 
//  - It contains a function 'eval<I>()', where I is a integer
//    parameter, 0, 1 or 2, which returns the value of the Ith element
//    of the vector expression. [For matrix expressions this is replaced
//    by 'eval<I,J>()'.]
// 
//  - The 'operator+' is overloaded to return an instance of this class
// 
//  - The 'operator= ' is overloaded to call the 'eval<I>()' function.
// 
//  - Two class template-parameters are needed to specify the operant types.
//    This way one can distinguish different operations, e.g. vector-vector 
//    multiplication from vector-double multiplication. 
// 
//  - There are separate expression templates Vector-valued and Matrix-valued 
//    expressions, such that one can distinguish e.g. Matrix-Vector from 
//    Matrix-Matrix multiplication.
// 
//  The general expression templates are defined as follows:
//    template <typename T, class A, operation B, class C>
//    Vector<T,A,B,C>; 
//  and
//    template <typename T, class A, operation B, class C>
//    Matrix<T,A,B,C>;
//  
//  The general template classes are not defined, only special instances for 
//  allowed operations B. 
// 
//  The default values for the templates arguments are 'Base' for 'A'
//  and 'C', and 'NoOp' for 'B'. The templates with these default
//  values, denoted as 'Vector<TT>' and 'Matrix<TT>' serve as the actual
//  Vector and Matrix class. They are partially specialized template
//  classes in the vecmat3 namespace, with the template parameter T
//  signifying the data type.  For convenience and backward
//  compatibility, Vector<DOUBLE> and Matrix<DOUBLE> and are
//  'typedef'ed as the default Vector and Matrix classes.  The
//  specializations of these default templates therefore contain the
//  actual elements and much of the basic functionality of vectors and
//  matrices.
// 
//  Note that the classes A and C can themselved be expression classes,
//  and so nested expressions such as (a+b)*(d+2*c) are perfectly
//  possible.  This recursiveness does make the notation somewhat
//  involved.
// 
//  For further information on the technique of expression templates
//  in c++, see:
//
// - ubiety.uwaterloo.ca/~tveldhui/papers/Expression-Templates/expre_eval_l.html
//
// - www.oonumerics.org/blitz
//
// - tvmet.sourceforge.net
//
