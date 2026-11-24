#define ASSERT( cond, msg )\
  if( .not. cond ) then; \
     write(*,"('> error at >:',1x,a,'::',i0,':')")__FILE__,__LINE__ ;\
     write(*,"('>',1x,*(g0))") msg ;\
     error stop 1 ;\
  end if

