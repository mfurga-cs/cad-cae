%% 
function spline2D()
  spline2Duniform()
return

function t=check_sanity(knot_vector,p)

  initial = knot_vector(1)
  kvsize = size(knot_vector,2)

  t = true
  counter = 1

  for i=1:p+1
    if (initial ~= knot_vector(i))
      t = false
      return
    end
  end


  for i=p+2:kvsize-p-1
    if (initial == knot_vector(i))
      counter = counter + 1
      if (counter > p)
        t = false  
        return
      end
    else
      initial = knot_vector(i)
      counter = 1
    end
  end

  initial = knot_vector(kvsize)

  for i=kvsize-p:kvsize
    if (initial ~= knot_vector(i))
      t = false
      return
    end
  end
  
  for i=1:kvsize-1
    if (knot_vector(i)>knot_vector(i+1))
      t = false
    end
  end

  return

end

function spline2Duniform()
  img = imread("aerial.jpg");
  img = rgb2gray(img);

  h = size(img, 1);
  w = size(img, 2);

  % [ -2 -1 0 1 2 3 ... w-3 w-2 w-1 w ]
  % w = 300
  % [ 0 0 0 1 2 3 ... 297 298 298 298]
  %
  % nrx = w + 2
  % nry = h + 2

  knot_vectorx = -2 : w;
  knot_vectorx(1:3) = 0;
  knot_vectorx(end - 1: end) = w - 2;

  knot_vectory = -2 : h;
  knot_vectory(1:3) = 0;
  knot_vectory(end - 2: end) = h - 2;

  coeff_vector = double(img) / 255;
  precision = 0.01

  % image(img);
  % return;


  %macros
  compute_nr_basis_functions = @(knot_vector,p) size(knot_vector, 2) - p - 1
  mesh   = @(a,c) [a:precision*(c-a):c]

  %splines in x
  px = compute_p(knot_vectorx)
  tx = check_sanity(knot_vectorx,px)
  nrx = compute_nr_basis_functions(knot_vectorx,px)

  x_begin = knot_vectorx(1)
  x_end = knot_vectorx(size(knot_vectorx,2))

  x=mesh(x_begin,x_end);
  %splines in y
  py = compute_p(knot_vectory)
  ty = check_sanity(knot_vectory,py)
  nry = compute_nr_basis_functions(knot_vectory,py)

  y_begin = knot_vectory(1)
  y_end = knot_vectory(size(knot_vectory,2))


  y=mesh(y_begin,y_end);

  %X and Y coordinates of points over the 2D mesh
  [X,Y]=meshgrid(x,y);

  hold on

  u = 0;

  for i=1:nrx
    vx=compute_spline(knot_vectorx,px,i,X);
    for j=1:nry
      vy=compute_spline(knot_vectory,py,j,Y);
      u = u + coeff_vector(h - j + 1, i) * vx .* vy;
    end
  end

  surf(X, Y, u);

  hold off
end

function p=compute_p(knot_vector)

  initial = knot_vector(1)
  kvsize = size(knot_vector,2)
  i=1

  while (i+1 < kvsize) && (initial == knot_vector(i+1))
    i=i+1
  end
  
  p = i-1
  
  return
end

function y=compute_spline(knot_vector,p,nr,x)
  
  fC= @(x,a,b) (x)/(b-a)-a/(b-a);
  fD= @(x,c,d) (1-x)/(d-c)+(d-1)/(d-c);
  
  
  a = knot_vector(nr);
  b = knot_vector(nr+p);
  c = knot_vector(nr+1);
  d = knot_vector(nr+p+1);

  if (p==0)
    y = 0 .* (x < a) + 1 .* (a <= x & x <= d) + 0 .* (x > d);
    return
  end
  
  lp = compute_spline(knot_vector,p-1,nr,x);
  rp = compute_spline(knot_vector,p-1,nr+1,x);
  
  if (a==b)
    y1 = 0 .* (x < a) + 1 .* (a <= x & x <= b) + 0 .* (x > b);
  else
    y1 = 0 .* (x < a) + fC(x,a,b) .* (a <= x & x <= b) + 0 .* (x > b);
  end
  
  if (c==d)
    y2 = 0 .* (x < c) + 1 .* (c < x & x <= d) + 0 .* (d < x);
  else
    y2 = 0 .* (x < c) + fD(x,c,d) .* (c < x & x <= d) + 0 .* (d < x);
  end
  
  y = lp .* y1 + rp .* y2;
  return
  
end

end
