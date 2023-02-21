function ress = mtimes(a,bb)

if a.adjoint,
    % Multicoil non-Cartesian k-space to Cartesian image domain
    % nufft for each coil and time point
    for tt=1:size(bb,4),
        for ch=1:size(bb,3),
            b = bb(:,:,ch,tt).*a.w(:,:,tt);
            res(:,:,ch) = reshape(nufft_adj(b(:),a.st{tt})/sqrt(prod(a.imSize2)),a.imSize(1),a.imSize(2));
        end
        ress(:,:,tt)=sum(res.*conj(a.b1),3)./sum(abs((squeeze(a.b1))).^2,3);
        clear res
    end
    ress=ress.*size(a.w,1)*pi/2/size(a.w,2);
else
    % Cartesian image to multicoil non-Cartesian k-space
%     [m, n] = size(unwrapped_phase);
    wrapped_phase = 2*atan(sin(unwrapped_phase)./(1 + cos(unwrapped_phase)));
    wrapped_phase(cos(unwrapped_phase) == -1)= pi;
%     for i = 1 : m
%         for j = 1 : n
%             if (cos(unwrapped_phase(i,j)) == -1)
%                 wrapped_phase(i,j) = pi;
%             end
%         end
%     end
end

