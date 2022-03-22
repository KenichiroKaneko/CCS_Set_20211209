load("cal/result_2033_LcurveKUP")
figure()
hold on 
plot(51:75, err_LCFS_list(51:75))
[m I] = min(err_LCFS_list)
plot(I, m, 'o')
