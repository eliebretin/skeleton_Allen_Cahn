function affiche_solution_3d2(x,U1_0,U)
clf;


v = U;
alpha(1)
 p = patch(isosurface(x,x,x,v,0.5));
 isonormals(x,x,x,v,p)
 set(p,'FaceColor','cyan','EdgeColor','none');

 
 
 v = U1_0;
  p = patch(isosurface(x,x,x,v,0.5));
  isonormals(x,x,x,v,p)
  set(p,'FaceColor','red','EdgeColor','none');
%alpha(0.9)
alpha(0.95)

daspect([1 1 1])
 view([1,0.5,1]); 
 %camlight headlight;
%camlight left;
%camlight(2,2)
 camlight('infinite')
lighting gouraud

 %view(3);
%view([150,40])
%view([100,20]) 
%view([90,20])
%view([200,60])
%view([95,10]) 
