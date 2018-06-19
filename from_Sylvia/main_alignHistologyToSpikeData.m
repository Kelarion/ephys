subject = 'SS087';
date = '2017-12-12';

folder = ['C:\DATA\Spikes\' subject '\' date '\histology'];

av = readNPY('C:\Users\Matteo A\Documents\MATLAB\AllenCCF\allenCCF\annotation_volume_10um_by_index.npy');
tv = readNPY('C:\Users\Matteo A\Documents\MATLAB\AllenCCF\allenCCF\template_volume_10um.npy');
st = loadStructureTree('structure_tree_safe_2017.csv');
f = allenAtlasBrowser(tv,av,st);

%% mark points and save them (before closing the atlas)!

tracedPoints = f.UserData.pointList(:, [3 2 1]);
if tracedPoints(1,3)>570
    tracedPoints(:,3) = 1140-tracedPoints(:,3);
end
save(fullfile(folder, sprintf('%s_%s_tracedPoints.mat', subject, date)), ...
    'tracedPoints')

[m,p,s] = best_fit_line(tracedPoints(:,1), tracedPoints(:,2), tracedPoints(:,3));
% ensure proper orientation: want 0 at the top of the brain and positive
% distance goes down into the brain (opposite of the probe, but works with
% scaling better)
if p(2)<0
    p = -p;
end
% determine "origin" at top of brain
ann = 10;
gotToCtx = false;
isoCtxId = num2str(st.id(strcmp(st.acronym, 'Isocortex')));
while ~(ann==1 && gotToCtx)
    m = m-p; % step 10um, backwards up the track
    ann = av(round(m(1)),round(m(2)),round(m(3)));
    if ~isempty(strfind(st.structure_id_path{ann}, isoCtxId))
        % if the track didn't get to cortex yet, keep looking, might be in
        % a gap between midbrain/thal/etc and cortex. Once you got to
        % isocortex once, you're good. 
        gotToCtx = true;
    end
end
fwireframe = plotBrainGrid([]);
% figure(fwireframe);
hold on; 
hp = plot3(tracedPoints(:,1), tracedPoints(:,3), tracedPoints(:,2), 'o');
plot3(m(1), m(3), m(2), 'k*')
t = 0:600;
plot3(m(1)+p(1)*t([1 end]), m(3)+p(3)*t([1 end]), m(2)+p(2)*t([1 end]), ...
    'Color', get(hp, 'Color'), 'LineWidth', 2);

xc = repmat([43;11;59;27], 150,1); % build coords without refs excluded, for visual clarity
yc = reshape(repmat(20*[0:299], 2,1), 600, 1);

[borders, fD, ~] = plotLabelsAsProbe(m, p, xc, yc-20, av, st, 1);
xlim([0 67])
ylim([0 6000])
set(gcf, 'Position', [1293 42 220 1074])