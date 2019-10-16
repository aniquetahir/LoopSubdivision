#include <fstream>
#include <vector>
#include <iostream>

#include "OBJFileReader.h"
#include "Solid.h"
#include "iterators.h"
#include <sstream>
#include <string>

#include <map>

using namespace std;
using namespace MeshLib;


template <class T>
vector<T*> removeVectorDuplicates(vector<T*> v){
    vector<T*> uniqueElementVector;
    for(int i=0; i<v.size(); i++){
        int firstIndex = i;
        for(int j=0; j<v.size(); j++){
            if(i!=j && v[i]==v[j]){
                firstIndex = min(firstIndex, j);
            }
        }
        if(firstIndex == i){
            uniqueElementVector.push_back(v[i]);
        }
    }

    return uniqueElementVector;
}

template <class T>
bool elementExists(vector<T*> v, T* element){
    for(auto ve:v){
        if(ve==element){
            return true;
        }
    }
    return false;
}

//map<Point*, HalfEdge**> getIntermediatePoints(vector<HalfEdge*> allHE){
//    map<Point*, HalfEdge**> intermediatePoints;
//    vector<HalfEdge*> usedHE;
//
//    for(auto he: allHE){
//        if(elementExists(usedHE, he)){
//            continue;
//        }
//        Point p1 = he->vertex()->point();
//        Point p2 = he->he_next()->vertex()->point();
//
//        Point* midPoint = new Point(
//                (p1[0]+p2[0])/2,
//                (p1[1]+p2[1])/2,
//                (p1[2]+p2[2])/2 );
//
//
//        HalfEdge** heArray = (HalfEdge**)malloc(sizeof(HalfEdge**) * 2);
//
//        heArray[0] = he;
//        heArray[1] = he->he_sym();
//
//        usedHE.push_back(he);
//        usedHE.push_back(he->he_sym());
//
//        intermediatePoints[midPoint] = heArray;
//    }
//
//    return intermediatePoints;
//}

string stringifyVertex(Vertex* v){
    std::stringstream ss;
    ss<<v->point()[0]<<","<<v->point()[1]<<","<<v->point()[2];
    return ss.str();
}

vector<Vertex*> getVertexNeighbors(Vertex* v){
    bool isBoundary = v->boundary();

    vector<Vertex*> nbrs;

    HalfEdge* origHe = v->halfedge();
    HalfEdge* he = v->halfedge();
    do{
        Vertex* nbr = he->source();
        if(!isBoundary){
            nbrs.push_back(nbr);
        }else if(he->he_sym() == NULL){
            nbrs.push_back(nbr);
        }
        if(he->he_next()->he_sym() == NULL && isBoundary){
            nbrs.push_back(he->he_next()->target());
            // go back to find other neighbor
            if(nbrs.size()<2){
                while (he->he_sym()!=NULL){
                    he=he->he_sym()->he_prev();
                }
                nbrs.push_back(he->source());
            }
            break;
        }

        he = he->he_next()->he_sym();
    }while(he!=origHe);

    return nbrs;
}


int main(int argc, char *argv[])
{
	// Read in the obj file
	Solid mesh;
	OBJFileReader of;
	std::ifstream in(argv[1]);
	of.readToSolid(&mesh, in);

    // Create new vertices at the center of each edge
    vector<HalfEdge*> allHalfEdges;
    vector<HalfEdge*> parsedHalfEdges;
    vector<Vertex*> oldPoints;
    vector<Vertex*> newPoints;

    //cout<< "x,y,z,color"<< std::endl;

	// Get all points
	SolidVertexIterator vIter(&mesh);
    for(;!vIter.end(); ++vIter){
        Vertex* v = *vIter;
        oldPoints.push_back(v);
        allHalfEdges.push_back(v->halfedge());
    }

    vector<Edge*> old_edges;
    SolidEdgeIterator edge_iter(&mesh);
    map<Edge*, vector<Vertex*>> edgeNeighbors;

    for(;!edge_iter.end(); ++edge_iter){
        Edge* e = *edge_iter;
        old_edges.push_back(e);

        Vertex* v1;
        Vertex* v2;
        e->get_vertices(v1, v2);

        vector<Vertex*> nPoints;
        nPoints.push_back(v1);
        nPoints.push_back(v2);

        if(!e->boundary()){
            nPoints.push_back(e->halfedge(0)->he_next()->target());
            nPoints.push_back(e->halfedge(1)->he_next()->target());
        }

        edgeNeighbors[e] = nPoints;
    }

    map<Vertex*, vector<Vertex*>> neighbors;
    vector<Edge*> tbs; // To Be Swapped

    // Populate Neighbours of new points, Create new faces/points
    for(auto edge: old_edges){
        Vertex* v1;
        Vertex* v2;
        edge->get_vertices(v1, v2);

        vector<Vertex*> nPoints = edgeNeighbors[edge];

        Point newPoint = (v1->point() + v2->point())/2;

        AddedElements elemAdded = mesh.edgeSplitGetNewEdges(edge);
        Vertex * newVertex = elemAdded.vertex;

        newVertex->point() = newPoint;

        newPoints.push_back(newVertex);
        neighbors[newVertex] = nPoints;

        for(long i=0;i<elemAdded.num_edges;i++){
            Edge* edge = elemAdded.edges[i];
            Vertex* eV1;
            Vertex* eV2;

            edge->get_vertices(eV1, eV2);
            if(!(elementExists(newPoints, eV1) && elementExists(newPoints, eV2))){
                tbs.push_back(edge);
            }
        }
    }


    // Now that we have all new vertices. We swap edges between a new vertex and an old one that is incorrect
    for(auto edge:tbs){
        mesh.edgeSwap(edge);
    }


    // Create the new positions according to the algorithm
    // For new/odd vertices
    for(auto pair:neighbors){
        Vertex* np = pair.first;
        vector<Vertex*> np_neighborhood = pair.second;

//        cout<<"Point: ("<<stringifyVertex(np)<<")"<<endl;
//        cout<<"Neighbors: ";
//        for (auto n: np_neighborhood){
//            cout<<"("<<stringifyVertex(n)<<"),";
//        }
//        cout<<endl<<"============"<<endl;

        if(np_neighborhood.size() == 2){ //
            np->point() = (np_neighborhood[0]->point() * 0.5) + (np_neighborhood[1]->point() * 0.5);
        }else{
            double three_by_eight = (double)3/(double)8;
            double one_by_eight = (double)1/(double)8;
            np->point() =
                    (np_neighborhood[0]->point() * three_by_eight)
                    + ((np_neighborhood[1]->point() * three_by_eight))
                    + (np_neighborhood[2]->point() * one_by_eight)
                    + (np_neighborhood[3]->point() * one_by_eight);
        }
    }

    // For old/even vertices
    for(auto op:oldPoints){
        vector<Vertex*> opNbrs = getVertexNeighbors(op);

//        cout<<"Point: ("<<stringifyVertex(op)<<")"<<endl;
//        cout<<"Neighbors: ";
//        for (auto n: opNbrs){
//            cout<<"("<<stringifyVertex(n)<<"),";
//        }


        if(op->boundary()){
            double one_by_eight = (double)1/(double)8;
            double three_by_four = (double)3/(double)4;
            op->point() = opNbrs[0]->point()*one_by_eight + opNbrs[1]->point()*one_by_eight + op->point()*three_by_four;
        }else{
            int num_nbrs = opNbrs.size();
            Point* p = &op->point();
            double beta = 0;
            if(num_nbrs==3){
                beta = (double)3/(double)16;
            }else{
                beta = (double)3/((double)8*(double)num_nbrs);
            }

            *p = *p*(1-num_nbrs*beta);
            for(auto nbr: opNbrs){
                *p = ((nbr->point()*beta)) + *p;
            }

        }
//        cout<<endl<<"New Position: ("<<stringifyVertex(op)<<")";
//
//        cout<<endl<<"============"<<endl;
    }


//    removeVectorDuplicates(allHalfEdges);
//    map<Point*, HalfEdge**> intermediatePoints = getIntermediatePoints(allHalfEdges);
//
//    map<HalfEdge*, HalfEdge**> halfedgeReplacements;
//
//    // TODO Calculate new vertex positions
//
//    // Create Half Edge replacements
//    for(auto he: allHalfEdges){
//
//        HalfEdge* newHE1 = new HalfEdge();
//        HalfEdge* newHE2 = new HalfEdge();
//
//        newHE1.
//
//    }

    // Create new faces

    // Remove previous edges OR create mesh from new half edges


//    for(auto he: allHalfEdges){
//        if(!elementExists(parsedHalfEdges, he)){
//            Point p1 = he->vertex()->point();
//            Point p2 = he->he_next()->vertex()->point();
//
//            Vertex * v = mesh.edgeSplit(he->edge());
//            v->point() = (p1+p2)/2;
//
//            newPoints.push_back(v);
//            parsedHalfEdges.push_back(he);
//            parsedHalfEdges.push_back(he->he_sym());
//        }
//    }
//
//
//
//    for(auto m_pair:intermediatePoints){
//        Point p = *(m_pair.first);
//        cout<<p[0] << "," << p[1] << "," << p[2] << "," << "blue" << std::endl;
//    }








    /******************* Put you subdivision processing here *********************/

	// Write out the resultant obj file
	int vObjID = 1;
	std::map<int, int> vidToObjID;

	std::ofstream os(argv[2]);

	SolidVertexIterator iter(&mesh);

	for(; !iter.end(); ++iter)
	{
		Vertex *v = *iter;
		Point p = v->point();
		os << "v " << p[0] << " " << p[1] << " " << p[2] << std::endl;
		vidToObjID[v->id()] = vObjID++;
	}
	os << "# " << (unsigned int)mesh.numVertices() << " vertices" << std::endl;

	float u = 0.0, v = 0.0;
	for(iter.reset(); !iter.end(); ++iter)
	{
		Vertex *vv = *iter;
		std::string key( "uv" );
		std::string s = Trait::getTraitValue (vv->string(), key );
		if( s.length() > 0 )
		{
			sscanf( s.c_str (), "%f %f", &u, &v );
		}
		os << "vt " << u << " " << v << std::endl;
	}
	os << "# " << (unsigned int)mesh.numVertices() << " texture coordinates" << std::endl;

	SolidFaceIterator fiter(&mesh);
	for(; !fiter.end(); ++fiter)
	{
		Face *f = *fiter;
		FaceVertexIterator viter(f);
		os << "f " ;
		for(; !viter.end(); ++viter)
		{
			Vertex *v = *viter;
			os << vidToObjID[v->id()] << "/" << vidToObjID[v->id()] << " ";
		}
		os << std::endl;
	}
	os.close();

	return 0;
}