29a30,41
>     /**
>        the face-index of the surface element maps into
>        this table.
>     */
>     Array<FaceDescriptor> facedecoding;
>     /// element -> face, element -> edge etc ...
>     class MeshTopology * topology;
> 
>     /// nodes identified by close points 
>     class AnisotropicClusters * clusters;
> 
>     
71,75d82
<     /**
<        the face-index of the surface element maps into
<        this table.
<     */
<     Array<FaceDescriptor> facedecoding;
102,103d108
<     /// element -> face, element -> edge etc ...
<     class MeshTopology * topology;
107,109d111
<     /// nodes identified by close points 
<     class AnisotropicClusters * clusters;
< 
469a472,473
>     ///
>     void Merge (Mesh* m, const int surfindex_offset = 0);
