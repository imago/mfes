142,143c142,160
< 
< 
---
>    
>    // Changes the bc number of a VOL mesh
>    DLL_HEADER Ng_Mesh * Ng_SetProperties(Ng_Mesh * mesh, int surfnr, int bcnr, int domin, int domout)
>    {
>       Mesh * m = (Mesh*) mesh;
>       int maxsteps = m->GetNSE();
>       SurfaceElementIndex sei;
>       for (sei = 0; sei < maxsteps; sei++)
>       {
>     	//  if ((*m)[sei].GetIndex()){
>     		  m->GetFaceDescriptor((*m)[sei].GetIndex ()).SetSurfNr (surfnr);
>     		  m->GetFaceDescriptor((*m)[sei].GetIndex ()).SetBCProperty (bcnr);
>     		  m->GetFaceDescriptor((*m)[sei].GetIndex ()).SetDomainIn (domin);
>     		  m->GetFaceDescriptor((*m)[sei].GetIndex ()).SetDomainOut (domout);
>     	//  }
>       	        	 
>       }
>       return ( (Ng_Mesh*)m );
>    }
190c207,345
<       return NG_ERROR;
---
> 
> 	      Ng_Result status = NG_OK;
> 	      Mesh * m = (Mesh*)mesh1;
> 	      Mesh * m2 = (Mesh*)mesh2;
> 
> 	      if(!m)
> 	      {
> 	         status = NG_ERROR;
> 	      }
> 
> 	      if(status == NG_OK)
> 	      {
> 	         const int num_pts = m->GetNP();
> 	         int surfindex_offset = 0;
> 
> 	         int i;
> 
> 	         int oldnp = m->GetNP();
> 	         int oldne = m->GetNSeg();
> 
> 	         for(SurfaceElementIndex si = 0; si < m->GetNSE(); si++)
> 	           for(int j=1; j<=(*m)[si].GetNP(); j++) (*m)[si].GeomInfoPi(j).trignum = -1;
> 
> 	         int max_surfnr = 0;
> 	         for (i = 1; i <= m->GetNFD(); i++)
> 	           max_surfnr = max2 (max_surfnr, m->GetFaceDescriptor(i).SurfNr());
> 	         max_surfnr++;
> 
> 	         if(max_surfnr < surfindex_offset) max_surfnr = surfindex_offset;
> 
> 	         
> 	         int maxsteps = m2->GetNSE();
> 	         SurfaceElementIndex sei;
> 	         for (sei = 0; sei < maxsteps; sei++)
> 	         {
> 	        	 int j;
> 	        	 int surfnr, bcp, domin, domout, nep = 3, faceind = 0;
> 	        	 
> 	        	 surfnr = m2->GetFaceDescriptor((*m2)[sei].GetIndex ()).SurfNr()+1;
> 	        	 bcp = m2->GetFaceDescriptor((*m2)[sei].GetIndex ()).BCProperty()+1;
> 	        	 domin = m2->GetFaceDescriptor((*m2)[sei].GetIndex ()).DomainIn()+1; 
> 	        	 domout = m2->GetFaceDescriptor((*m2)[sei].GetIndex ()).DomainOut()+1; 
> 
> 	        	 for (j = 1; j <= m->facedecoding.Size(); j++)
> 	        		 if (m->GetFaceDescriptor(j).SurfNr() == surfnr &&
> 	        				 m->GetFaceDescriptor(j).BCProperty() == bcp &&
> 	        				 m->GetFaceDescriptor(j).DomainIn() == domin &&
> 	        				 m->GetFaceDescriptor(j).DomainOut() == domout)
> 	                             faceind = j;
> 
> 	              	  if (!faceind)
> 	              	  {
> 	              		  faceind = m->AddFaceDescriptor (FaceDescriptor(surfnr, domin, domout, 0));
> 	              		  if(m->GetDimension() == 2) bcp++;
> 	              		  m->GetFaceDescriptor(faceind).SetBCProperty (bcp);
> 	              	  }
> 
> 	                  Element2d tri(nep);
> 	                  tri.SetIndex(faceind);
> 
> 	                  Element2d sel = (*m2)[sei];
> 	                  
> 	                  for (j = 1; j <= nep; j++)
> 	                  {
> 	                	  tri.PNum(j) = sel.PNum(j) + oldnp;
> 	                  }
> 
> 
> 	                  for (j = 1; j <= nep; j++)
> 	                  {
> 	                	  tri.GeomInfoPi(j).trignum = -1;
> 	                  }
> 
> 	                  m->AddSurfaceElement (tri);
> 	         }
> 	        
> 	          for (ElementIndex ei = 0; ei < m2->GetNE(); ei++)
> 	          {
> 	        	  Element el;
> 	        	  
> 	        	  int hi = (*m2)[ei].GetIndex() + 1;
> 	        	  if (hi == 0) hi = 1;
> 	        	  el.SetIndex(hi);
> 	        	  
> 	        	  int nep = (*m2)[ei].GetNP();
> 	        	  el.SetNP(nep);
> 
> 	        	  for (int j = 0; j < nep; j++)
> 	        		  el[j] = (*m2)[ei][j] + oldnp;
> 
> 	        	  m->AddVolumeElement (el);
> 	         }
> 
> 	         
> 	         PointIndex pi;
> 	         maxsteps = m2->GetNP();
> 	         for (pi = PointIndex::BASE; pi < maxsteps+PointIndex::BASE; pi++)
> 	         {
> 	        	 Point3d p;
> 	        	 p.X() = (*m2)[pi](0)/1.0;
> 	        	 p.Y() = (*m2)[pi](1)/1.0;
> 	        	 p.Z() = (*m2)[pi](2)/1.0;
> 	             m->AddPoint (p);
> 	        	 
> 	         }  
> 	         
> 	         for (unsigned int i = 1; i <= m2->GetNSeg(); i++){
> 	        	 Segment & seg = m2->LineSegment (i);
> 	        	 
> 	             seg.surfnr1--;
>                  seg.surfnr2--;
>                  if(seg.surfnr1 >= 0)  seg.surfnr1 = seg.surfnr1 + max_surfnr;
>                  if(seg.surfnr2 >= 0)  seg.surfnr2 = seg.surfnr2 + max_surfnr;
>                  seg[0] = seg[0] +oldnp;
>                  seg[1] = seg[1] +oldnp;
>                  seg.edgenr = seg.edgenr + oldne;
>                  seg.epgeominfo[1].edgenr = seg.epgeominfo[1].edgenr + oldne;
> 	             m->AddSegment (seg);
> 	         }
> 	         
> 
> 	         m->CalcSurfacesOfNode ();
> 
> 	         m->topology -> Update();
> 	         m->clusters -> Update();
> 
> 	         m->SetNextMajorTimeStamp();
> 
> 	         if(m->GetNP() > num_pts)
> 	         {
> 	            status = NG_OK;
> 	         }
> 	         else
> 	         {
> 	            status = NG_ERROR;
> 	         }
> 	      }
> 	      return status;
> 
