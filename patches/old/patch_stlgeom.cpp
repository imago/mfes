43a44
>   ref = NULL;
63a65
>   delete ref;
101c103,106
<   return RefinementSTLGeometry (*this);
---
>   delete ref;
>   ref = new RefinementSTLGeometry(*this);
>   return *ref;
> 
2364,2368c2369,2370
<   int i;
<   for (i = 1; i <= NOExternalEdges(); i++)
<     {
<       AddEdge(GetExternalEdge(i).i1,GetExternalEdge(i).i2);
<     }
---
>   for (int i = 1; i <= NOExternalEdges(); i++)
>     AddEdge(GetExternalEdge(i).i1,GetExternalEdge(i).i2);
2496,2497c2498
<   int i;
<   for (i = 1; i <= edgedata->Size(); i++)
---
>   for (int i = 1; i <= edgedata->Size(); i++)
2523c2524
< 	  for (i = 1; i <= edgedata->Size(); i++)
---
> 	  for (int i = 1; i <= edgedata->Size(); i++)
2544c2545
<   for (i = 1; i <= edgedata->Size(); i++)
---
>   for (int i = 1; i <= edgedata->Size(); i++)
2628,2629c2629
<   int i;
<   for (i = 1; i <= GetNE(); i++)
---
>   for (int i = 1; i <= GetNE(); i++)
2651,2652c2651
<   int i,j;
<   for (i = 1; i <= GetNOFaces(); i++)
---
>   for (int i = 1; i <= GetNOFaces(); i++)
2658c2657
<   for (i = 1; i <= GetNT(); i++)
---
>   for (int i = 1; i <= GetNT(); i++)
2662c2661
<       for (j = 1; j <= 3; j++)
---
>       for (int j = 1; j <= 3; j++)
2668c2667
<   for (i = 1; i <= GetNOFaces(); i++)
---
>   for (int i = 1; i <= GetNOFaces(); i++)
2674,2675c2673,2674
<   int k, ap1, ap2;
<   for (i = 1; i <= GetNOFaces(); i++)
---
>   int ap1, ap2;
>   for (int i = 1; i <= GetNOFaces(); i++)
2680c2679,2682
< 	for (j = 1; j <= c.GetNChartT(); j++)
---
>         // bool foundone = false;
>         int longest_ap1, longest_ap2 = -1;
>         double maxlen = -1;
> 	for (int j = 1; j <= c.GetNChartT(); j++)
2683c2685
< 	    for (k = 1; k <= 3; k++)
---
> 	    for (int k = 1; k <= 3; k++)
2689c2691,2698
< 		    AddEdge(ap1,ap2);
---
>                     // AddEdge(ap1,ap2);
>                     double len = Dist(GetPoint(ap1), GetPoint(ap2));
>                     if (len > maxlen) 
>                       {
>                         maxlen = len;
>                         longest_ap1 = ap1;
>                         longest_ap2 = ap2;
>                       }
2693a2703,2704
>         if (maxlen > 0)
>           AddEdge(longest_ap1,longest_ap2);
