415c415
<             outfile << " " << GetFaceDescriptor((*this)[sei].GetIndex ()).SurfNr()+1;
---
>             outfile << " " << GetFaceDescriptor((*this)[sei].GetIndex ()).SurfNr();
1244,1245c1244,1248
<                 if(domin > 0) domin += oldnd;
<                 if(domout > 0) domout += oldnd;
---
> 		 	 // if(domin > 0)  domin += oldnd;
> 				domin += oldnd;
>        	     // if(domout > 0)  domout += oldnd;
> 				domout += oldnd;
> 
