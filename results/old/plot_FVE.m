clc;
figure(250);
stem(max(m1(kv_locs))-m1(kv_locs),abs(b1(kv_locs)))
ylabel('|B1|');
xlabel('K_v (s/cm)');

vfov=(1./diff(m1(kv_locs)))/42.58e6;
disp('velocity FOV:')
disp(vfov);

delta_kv=(diff(m1(kv_locs)));
disp('Delta K_v:')
disp(delta_kv);
