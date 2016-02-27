using UnityEngine;
using System.Collections;
using System.Collections.Generic;
using TriangleNet;
using TriangleNet.Geometry;
using TriangleNet.Tools;

public class TriangleConverter {

    public Vector3[] Vertices;
    public int[] Triangles;
    /*public Vector2[] uv;*/

    public void SetInputGeometry(InputGeometry data)
    {

    }

    public void SetMesh(TriangleNet.Mesh mesh)
    {
        int n = mesh.Vertices.Count;
        int i = 0;

        // Linear numbering of mesh
        mesh.Renumber();

        // Copy points
        this.Vertices = new Vector3[n];
        foreach(var pt in mesh.Vertices)
        {
            this.Vertices[i] = new Vector3((float)pt.X, (float)pt.Y, 0);
            i++;
        }

        // Copy Triangles
        var triangles = new List<int>(3 * mesh.Triangles.Count);
        foreach (var tri in mesh.Triangles)
        {
            triangles.Add(tri.P0);
            triangles.Add(tri.P1);
            triangles.Add(tri.P2);

            /*if (this.NumberOfRegions > 0)
            {
                this.TrianglePartition[i++] = tri.Region;
            }*/
        }
        this.Triangles = triangles.ToArray();

    }
}
