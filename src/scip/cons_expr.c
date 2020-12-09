#ifdef SCIP_DISABLED_CODE
/** compares nonlinear handler by enforcement priority
 *
 * if handlers have same enforcement priority, then compare by detection priority, then by name
 */
static
int nlhdlrEnfoCmp(
   void*                 hdlr1,              /**< first handler */
   void*                 hdlr2               /**< second handler */
)
{
   SCIP_NLHDLR* h1;
   SCIP_NLHDLR* h2;

   assert(hdlr1 != NULL);
   assert(hdlr2 != NULL);

   h1 = (SCIP_NLHDLR*)hdlr1;
   h2 = (SCIP_NLHDLR*)hdlr2;

   if( h1->enfopriority != h2->enfopriority )
      return (int)(h1->enfopriority - h2->enfopriority);

   if( h1->detectpriority != h2->detectpriority )
      return (int)(h1->detectpriority - h2->detectpriority);

   return strcmp(h1->name, h2->name);
}
#endif
