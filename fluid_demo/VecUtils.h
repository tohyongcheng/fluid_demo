#ifndef VEC_UTILS_H
#define VEC_UTILS_H

#include "vecmath/include/vecmath.h"

class VecUtils
{
public:
    
    static Vector3f min( const Vector3f& b, const Vector3f& c )
    {
        Vector3f out;
        
        for( int i = 0; i < 3; ++i )
        {
            out[ i ] = ( b[i] < c[i] ) ? b[i] : c[i];
        }
        
        return out;
    }
    
    static Vector3f max( const Vector3f& b, const Vector3f& c )
    {
        Vector3f out;
        
        for( int i = 0; i < 3; ++i )
        {
            out[ i ] = ( b[i] > c[i] ) ? b[i] : c[i];
        }
        
        return out;
    }
    
    static Vector3f clamp( const Vector3f& data, float low = 0, float high = 1 )
    {
        Vector3f out = data;
        for( int i = 0; i < 3; ++i )
        {
            if( out[ i ] < low )
            {
                out[ i ] = low;
            }
            if( out[ i ] > high )
            {
                out[ i ] = high;
            }
        }
        
        return out;
    }
    
    // transforms a 3D point using a matrix, returning a 3D point
    static Vector3f transformPoint( const Matrix4f& mat, const Vector3f& point )
    {
        return ( mat * Vector4f( point, 1 ) ).xyz();
    }
    
    // transform a 3D directino using a matrix, returning a direction
    // This function *does not* take the inverse tranpose for you.
    static Vector3f transformDirection( const Matrix4f& mat, const Vector3f& dir )
    {
        return ( mat * Vector4f( dir, 0 ) ).xyz();
    }
    
    static Vector3f cramer( const Matrix3f A, const Vector3f B) {
        Vector3f result;
        float det_a = A.determinant();
        assert(det_a!=0);
        result[0] = Matrix3f(B, A.getCol(1), A.getCol(2)).determinant()/det_a;
        result[1] = Matrix3f(A.getCol(0), B, A.getCol(2)).determinant()/det_a;
        result[2] = Matrix3f(A.getCol(0), A.getCol(1), B).determinant()/det_a;
        return result;
    }
    
};

#endif // VEC_UTILS_H
