256a257,270
> /*! \brief Set element options
> 
>     Changes the element options of a Netgen mesh in VOL format . 
> 	(Sakalli)
> 	
>     \param mesh       Name of the Mesh structure already existent in memory
>     \param surfnr	  Set surface number
>     \param bcnr 	  Set boundary number
>     \param domin 	  Set domain in
>     \param domout 	  Set domain out
>     \return Ng_Result Status of the boundary operation
> */
> DLL_HEADER Ng_Mesh * Ng_SetProperties(Ng_Mesh * m, int surfnr, int bcnr, int domin, int domout);
> 
274c288
<     (NOTE: FUNCTION STILL WORK IN PROGRESS!!!)
---
>     (NOTE: First version implemented, Sakalli)
