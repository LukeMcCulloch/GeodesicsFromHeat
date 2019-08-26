#include "DiscreteExteriorCalculus.h"

namespace DDG
{


/*

      Build a sparse matrix encoding the Hodge operator on 0-forms for this mesh.
      Returns a sparse, diagonal matrix corresponding to vertices.

      The discrete hodge star is a diagonal matrix where each entry is
      the (area of the dual element) / (area of the primal element). You will probably
      want to make use of the Vertex.circumcentricDualArea property you just defined.

      TLM as seen in notes:
      By convention, the area of a vertex is 1.0.
*/
   template <class T>
   void HodgeStar0Form<T> :: build( const Mesh& mesh,
                                    SparseMatrix<T>& star0 )
   // builds a diagonal matrix mapping primal discrete 0-forms
   // to dual discrete 2-forms
   {
      int nV = mesh.vertices.size();

      star0 = SparseMatrix<T>( nV, nV );

      for( VertexCIter v  = mesh.vertices.begin();
                       v != mesh.vertices.end();
                       v ++ )
      {
         int i = v->index;
         star0( i, i ) = v->area();
      }
   }

/*

      Build a sparse matrix encoding the Hodge operator on 1-forms for this mesh.
      Returns a sparse, diagonal matrix corresponding to edges.

      The discrete hodge star is a diagonal matrix where each entry is
      the (area of the dual element) / (area of the primal element). The solution
      to exercise 26 from the previous homework will be useful here.

      TLM: cotan formula again.  see ddg notes page 89
      see also source slide 56:
         http://brickisland.net/DDGFall2017/wp-content/uploads/2017/09/
                              CMU_DDG_Fall2017_06_DiscreteExteriorCalculus.pdf

      Note that for some geometries, some entries of hodge1 operator may be exactly 0.
      This can create a problem when we go to invert the matrix. To numerically sidestep
      this issue, you probably want to add a small value (like 10^-8) to these diagonal 
      elements to ensure all are nonzero without significantly changing the result.

*/
   template <class T>
   void HodgeStar1Form<T> :: build( const Mesh& mesh,
                                    SparseMatrix<T>& star1 )
   // builds a diagonal matrix mapping primal discrete 1-forms
   // to dual discrete 1-forms
   {
      int nE = mesh.edges.size();

      star1 = SparseMatrix<T>( nE, nE );

      for( EdgeCIter e  = mesh.edges.begin();
                     e != mesh.edges.end();
                     e ++ )
      {
         // get the cotangents of the two angles opposite this edge
         double cotAlpha = e->he->cotan();
         double cotBeta  = e->he->flip->cotan();

         int i = e->index;
         star1( i, i ) = ( cotAlpha + cotBeta ) / 2.;
      }
   }

   /*

         Build a sparse matrix encoding the Hodge operator on 2-forms for this mesh
         Returns a sparse, diagonal matrix corresponding to faces.

         The discrete hodge star is a diagonal matrix where each entry is
         the (area of the dual element) / (area of the primal element).


         hint hint!, vertex is => (dual) vertex:
         By convention, the area of a vertex is 1.0.


         TLM: see also source slide 57:
            http://brickisland.net/DDGFall2017/wp-content/uploads/2017/09/
                              CMU_DDG_Fall2017_06_DiscreteExteriorCalculus.pdf
   */
   template <class T>
   void HodgeStar2Form<T> :: build( const Mesh& mesh,
                                    SparseMatrix<T>& star2 )
   // builds a diagonal matrix mapping primal discrete 2-forms
   // to dual discrete 2-forms
   {
      int nF = mesh.faces.size();

      star2 = SparseMatrix<T>( nF, nF );

      for( FaceCIter f  = mesh.faces.begin();
                     f != mesh.faces.end();
                     f ++ )
      {
         int i = f->index;
         star2( i, i ) = 1. / f->area();
      }
   }

   /*

         Build a sparse matrix encoding the exterior derivative on 0-forms.

         See section 3.6 of the course notes for an explanation of DEC.
        
         the exterior derivative of 0-form = [edge x vertex] adjacency matrix
   */
   template< class T >
   void ExteriorDerivative0Form<T> :: build( const Mesh& mesh,
                                             SparseMatrix<T>& d0 )
   {
      int nV = mesh.vertices.size();
      int nE = mesh.edges.size();

      d0 = SparseMatrix<T>( nE, nV );

      for( EdgeCIter e  = mesh.edges.begin();
                     e != mesh.edges.end();
                     e ++ )
      {
         // the row index is the index of the edge
         int r = e->index;
         
         // the column indices are the indices of the two
         // edge vertices -- orientation is determined by
         // the orientation of the edge's first half edge
         int ci = e->he->vertex->index;
         int cj = e->he->flip->vertex->index;

         d0( r, ci ) = -1.;
         d0( r, cj ) =  1.;
      }
   }


   /*

         Build a sparse matrix encoding the exterior derivative on 1-forms.

         See section 3.6 of the course notes for an explanation of DEC.
        
         the exterior derivative of 1-form = [face x edge] adjacency matrix
   */
   template< class T >
   void ExteriorDerivative1Form<T> :: build( const Mesh& mesh,
                                             SparseMatrix<T>& d1 )
   {
      int nE = mesh.edges.size();
      int nF = mesh.faces.size();

      d1 = SparseMatrix<T>( nF, nE );

      // visit each face
      for( FaceCIter f  = mesh.faces.begin();
                     f != mesh.faces.end();
                     f ++ )
      {
         // the row index is the index of the face
         int r = f->index;

         // visit all edges of this face
         HalfEdgeCIter he = f->he;
         do
         {
            // the column index is the index of the current edge
            int c = he->edge->index;

            // relative orientation is determined by checking if
            // the current half edge is the first half edge of its
            // corresponding edge
            double s = ( he->edge->he == he ? 1. : -1. );

            // set the entry for this edge
            d1( r, c ) = s;
            
            he = he->next;
         }
         while( he != f->he );
      }
   }
}

