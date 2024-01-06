//
// `#1#1#1, I am here.
//
   double Gradval[NbrPoints];
   u.GetGradients	(	ElementTransformation & 	tr,
const IntegrationRule & 	ir,
DenseMatrix & 	grad 
)	
   ElementTransformation *et;
   et = new(ElementTransformation);
   // get the value of each point one by one.
   for(i=0; i< NbrPoints; i++) {
   
      
      et.SetIntPoint	(&ips[i]);

     // u.GetGradient	(	&et, Vector & 	grad )	

   }