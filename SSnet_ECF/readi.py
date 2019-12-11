def readin(f):
  f1=open(f,`r`)
  fa=readlines()
  f1.close()
  fa0=fa[0]
  fd=fa0.split(`\t`)
  plist=[]
  for a in fa[1:]:
    b=a.split(`\t`)
    if b!=[]:
      plist.append(b)
  return([fd,plist])