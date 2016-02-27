using UnityEngine;
using System.Collections;
using TriangleNet.Geometry;
using TriangleNet.IO;

public class TriangleImporter : MonoBehaviour {

    TriangleNet.Mesh TriMesh;
	// Use this for initialization
	void Start () {

        TriMesh = new TriangleNet.Mesh();
        /*_mesh.Triangulate("Assets/Plugins/Data/superior.poly");*/

        var geometry = FileReader.ReadPolyFile("Assets/Plugins/Data/superior.poly");
        TriMesh.behavior.Quality = true;
        TriMesh.behavior.MinAngle = 30;
        TriMesh.Triangulate(geometry);

        Debug.Log("Size : " + TriMesh.Vertices.Count);

        int i = 0;
        foreach (var pt in TriMesh.Vertices)
        {
            Debug.Log(i++ + " : " + pt.X + "," + pt.Y);
        }

        foreach(var tri in TriMesh.Triangles)
        {
            /*Debug.Log(tri.)*/
        }

        FileWriter.Write(TriMesh, "Assets/Plugins/Data/superior_mesh.ele"); 
	}
	
}
