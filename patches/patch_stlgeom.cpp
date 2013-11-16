2364,2368c2364,2365
<   int i;
<   for (i = 1; i <= NOExternalEdges(); i++)
<     {
<       AddEdge(GetExternalEdge(i).i1,GetExternalEdge(i).i2);
<     }
---
>   for (int i = 1; i <= NOExternalEdges(); i++)
>     AddEdge(GetExternalEdge(i).i1,GetExternalEdge(i).i2);
2496,2497c2493
<   int i;
<   for (i = 1; i <= edgedata->Size(); i++)
---
>   for (int i = 1; i <= edgedata->Size(); i++)
2523c2519
< 	  for (i = 1; i <= edgedata->Size(); i++)
---
> 	  for (int i = 1; i <= edgedata->Size(); i++)
2544c2540
<   for (i = 1; i <= edgedata->Size(); i++)
---
>   for (int i = 1; i <= edgedata->Size(); i++)
2628,2629c2624
<   int i;
<   for (i = 1; i <= GetNE(); i++)
---
>   for (int i = 1; i <= GetNE(); i++)
2651,2652c2646
<   int i,j;
<   for (i = 1; i <= GetNOFaces(); i++)
---
>   for (int i = 1; i <= GetNOFaces(); i++)
2658c2652
<   for (i = 1; i <= GetNT(); i++)
---
>   for (int i = 1; i <= GetNT(); i++)
2662c2656
<       for (j = 1; j <= 3; j++)
---
>       for (int j = 1; j <= 3; j++)
2668c2662
<   for (i = 1; i <= GetNOFaces(); i++)
---
>   for (int i = 1; i <= GetNOFaces(); i++)
2674,2675c2668,2669
<   int k, ap1, ap2;
<   for (i = 1; i <= GetNOFaces(); i++)
---
>   int ap1, ap2;
>   for (int i = 1; i <= GetNOFaces(); i++)
2680c2674,2677
< 	for (j = 1; j <= c.GetNChartT(); j++)
---
>         // bool foundone = false;
>         int longest_ap1, longest_ap2 = -1;
>         double maxlen = -1;
> 	for (int j = 1; j <= c.GetNChartT(); j++)
2683c2680
< 	    for (k = 1; k <= 3; k++)
---
> 	    for (int k = 1; k <= 3; k++)
2689c2686,2693
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
2693a2698,2699
>         if (maxlen > 0)
>           AddEdge(longest_ap1,longest_ap2);
